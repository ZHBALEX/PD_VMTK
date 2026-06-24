from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import vtk

from .config import CenterlineConfig


@dataclass(slots=True)
class MeasurementResult:
    distances: list[float] = field(default_factory=list)
    points: list[np.ndarray] = field(default_factory=list)
    cross_section_areas: list[float] = field(default_factory=list)
    cross_section_perimeters: list[float] = field(default_factory=list)
    surface_areas: list[float] = field(default_factory=list)
    cross_sections: list[vtk.vtkPolyData] = field(default_factory=list)


def load_stl(path: str | Path) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(path))
    reader.Update()
    surface = reader.GetOutput()
    if surface is None or surface.GetNumberOfPoints() == 0:
        raise ValueError(f"Could not load STL surface: {path}")
    return surface


def clean_surface(surface: vtk.vtkPolyData) -> vtk.vtkPolyData:
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(surface)
    cleaner.Update()
    return cleaner.GetOutput()


def generate_centerline(surface, source_point, target_point, config: CenterlineConfig) -> vtk.vtkPolyData:
    method = config.centerline_method.lower()
    if method not in {"section_centroid", "centroid", "vmtk", "vtk_surface"}:
        raise ValueError(f"Unknown centerline method: {config.centerline_method}")

    if method == "section_centroid":
        return generate_section_centroid_centerline(surface, source_point, target_point, config)
    if method == "centroid":
        return generate_centroid_centerline(surface, source_point, target_point, config)
    if method == "vmtk":
        return generate_vmtk_centerline(surface, source_point, target_point, config)
    return generate_surface_path_centerline(surface, source_point, target_point, config)


def generate_vmtk_centerline(surface, source_point, target_point, config: CenterlineConfig) -> vtk.vtkPolyData:
    try:
        from vmtk import vmtkscripts
    except ImportError as exc:
        raise RuntimeError(
            "VMTK is required for proper duct/vessel centerline extraction. "
            "The VTK surface-path fallback is not a true centerline; choose vtk_surface only for debugging."
        ) from exc

    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = config.seed_selector_name
    centerlines.SourcePoints = [float(v) for v in source_point]
    centerlines.TargetPoints = [float(v) for v in target_point]
    centerlines.AppendEndPoints = config.append_end_points
    centerlines.Resampling = config.resampling
    centerlines.ResamplingStepLength = config.resampling_step_length
    centerlines.Execute()

    polydata = centerlines.Centerlines
    if polydata is None or polydata.GetNumberOfPoints() == 0:
        raise ValueError("Centerline generation failed. Check the model and selected points.")
    return orient_centerline(smooth_centerline(polydata, config.spline_filter_length), source_point)


def generate_centroid_centerline(surface, source_point, target_point, config: CenterlineConfig) -> vtk.vtkPolyData:
    cleaned = clean_surface(surface)
    surface_points = vtk_points_to_array(cleaned.GetPoints())
    if len(surface_points) < 3:
        raise ValueError("The surface does not contain enough points for centroid centerline extraction.")

    source = np.array(source_point, dtype=float)
    target = np.array(target_point, dtype=float)
    axis = target - source
    length = float(np.linalg.norm(axis))
    if length <= 0:
        raise ValueError("Source and target points are identical.")
    axis /= length

    progress = (surface_points - source) @ axis
    margin = max(length * 0.02, config.resampling_step_length)
    mask = (progress >= -margin) & (progress <= length + margin)
    points = surface_points[mask]
    progress = progress[mask]
    if len(points) < 3:
        raise ValueError("Could not find enough surface points between the selected endpoints.")
    if len(points) > 60000:
        sample_ids = np.linspace(0, len(points) - 1, 60000, dtype=int)
        points = points[sample_ids]
        progress = progress[sample_ids]

    step = effective_resampling_step(config, length / 500.0)
    smoothing_scale = max(config.spline_filter_length, step)
    n_samples = int(np.clip(length / max(step * 8.0, length / 96.0), 32, 120))
    centers = initial_centroid_samples(points, progress, length, n_samples)
    if len(centers) < 3:
        raise ValueError("Centroid centerline failed; try selecting endpoints farther apart.")

    center_points = centers
    smoothness = float(np.clip(0.45 + 0.45 * smoothing_scale / (smoothing_scale + step * 8.0), 0.45, 0.9))
    smooth_passes = int(np.clip(round(smoothing_scale / step), 2, 24))
    smooth_window = int(np.clip(round(smoothing_scale / step) | 1, 5, 21))
    center_points = trim_degenerate_polyline_ends(center_points, step)
    center_points = refine_principal_curve(points, center_points, passes=12 + smooth_passes, smoothness=smoothness)
    center_points = smooth_point_sequence(center_points, passes=smooth_passes, window=smooth_window)
    polydata = polyline_from_points(center_points)
    finished = finish_centerline(polydata, config, source_point)
    return enforce_monotone_selected_path(finished, source_point, target_point, step * 2.0)


def generate_section_centroid_centerline(surface, source_point, target_point, config: CenterlineConfig) -> vtk.vtkPolyData:
    cleaned = clean_surface(surface)
    guide = generate_centroid_centerline(cleaned, source_point, target_point, config)
    step = effective_resampling_step(config)
    smoothing_scale = max(config.spline_filter_length, step)
    guide_points = vtk_points_to_array(guide.GetPoints())
    if len(guide_points) < 3:
        raise ValueError("Guide path is too short for section-centroid centerline extraction.")

    guide_length = polyline_length(guide)
    sample_step = max(step * 6.0, guide_length / 140.0, 1e-6)
    guide = resample_polyline(guide, sample_step)
    guide_points = vtk_points_to_array(guide.GetPoints())

    center_points = guide_points
    for _ in range(4):
        center_points = refine_by_section_centers(cleaned, center_points, config, step)
        smooth_passes = int(np.clip(round(smoothing_scale / step), 2, 18))
        smooth_window = int(np.clip(round(smoothing_scale / step) | 1, 5, 17))
        center_points = smooth_point_sequence(center_points, passes=smooth_passes, window=smooth_window)
        center_points = trim_degenerate_polyline_ends(center_points, step)
        if len(center_points) < 3:
            return guide

    polydata = polyline_from_points(center_points)
    refined = finish_centerline(polydata, config, source_point)
    if refined.GetNumberOfPoints() < 2 or polyline_length(refined) < 0.6 * polyline_length(guide):
        return guide
    return enforce_monotone_selected_path(refined, source_point, target_point, step * 2.0)


def refine_by_section_centers(surface, guide_points: np.ndarray, config: CenterlineConfig, step: float) -> np.ndarray:
    centers = []
    last_center = None
    guide_list = [np.array(p, dtype=float) for p in guide_points]
    for i, guide_point in enumerate(guide_list):
        tangent = tangent_at(guide_list, i)
        center, area = section_center_and_area(surface, guide_point, tangent, config, last_center)
        if center is None:
            if last_center is None:
                continue
            center = last_center
        centers.append(center)
        last_center = center
    if len(centers) < 3:
        return guide_points
    centers = np.array(centers, dtype=float)
    return trim_degenerate_polyline_ends(centers, step)


def trim_degenerate_polyline_ends(points: np.ndarray, step: float) -> np.ndarray:
    if len(points) <= 4:
        return points
    cleaned = points.copy()
    min_move = max(step * 0.25, 1e-8)
    start = 0
    while start + 2 < len(cleaned) and np.linalg.norm(cleaned[start + 1] - cleaned[start]) < min_move:
        start += 1
    end = len(cleaned)
    while end - 3 >= start and np.linalg.norm(cleaned[end - 1] - cleaned[end - 2]) < min_move:
        end -= 1
    if end - start < 4:
        return points
    return cleaned[start:end]


def section_center(surface, guide_point: np.ndarray, tangent: np.ndarray, config: CenterlineConfig, last_center) -> np.ndarray | None:
    center, _area = section_center_and_area(surface, guide_point, tangent, config, last_center)
    return center


def section_center_and_area(surface, guide_point: np.ndarray, tangent: np.ndarray, config: CenterlineConfig, last_center) -> tuple[np.ndarray | None, float]:
    plane = vtk.vtkPlane()
    plane.SetOrigin(guide_point)
    plane.SetNormal(tangent)

    sphere = vtk.vtkSphere()
    sphere.SetCenter(guide_point)
    sphere.SetRadius(effective_section_radius(surface, config))

    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(surface)
    clipper.SetClipFunction(sphere)
    if config.sphere_inside_out:
        clipper.InsideOutOn()
    else:
        clipper.InsideOutOff()
    clipper.Update()

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(clipper.GetOutput())
    cutter.Update()

    stripper = vtk.vtkStripper()
    stripper.SetInputData(cutter.GetOutput())
    stripper.JoinContiguousSegmentsOn()
    stripper.Update()
    section = stripper.GetOutput()
    if section.GetNumberOfPoints() < 3 or section.GetNumberOfCells() == 0:
        return None, 0.0

    component, area = choose_section_component(section, guide_point, tangent, last_center)
    if component is None or len(component) < 3 or area <= 0.0:
        return None, 0.0
    return polygon_centroid_3d(np.array(component, dtype=float), tangent), area


def choose_section_component(section: vtk.vtkPolyData, guide_point: np.ndarray, tangent: np.ndarray, last_center) -> tuple[list[tuple[float, float, float]] | None, float]:
    ids = vtk.vtkIdList()
    best_score = float("inf")
    best_points = None
    best_area = 0.0
    reference = np.array(guide_point if last_center is None else last_center, dtype=float)
    for cell_id in range(section.GetNumberOfCells()):
        section.GetCellPoints(cell_id, ids)
        if ids.GetNumberOfIds() < 3:
            continue
        pts = [section.GetPoint(ids.GetId(i)) for i in range(ids.GetNumberOfIds())]
        arr = np.array(pts, dtype=float)
        area = abs(polygon_area_2d(project_to_plane(arr, tangent)))
        if area <= 1e-10:
            continue
        center = polygon_centroid_3d(arr, tangent)
        nearest_to_guide = np.min(np.linalg.norm(arr - guide_point, axis=1))
        continuity = np.linalg.norm(center - reference)
        radius = np.sqrt(area / np.pi)
        score = 0.75 * nearest_to_guide + 1.2 * continuity - 0.15 * radius
        if score < best_score:
            best_score = score
            best_points = pts
            best_area = area
    return best_points, best_area


def initial_centroid_samples(points: np.ndarray, progress: np.ndarray, length: float, n_samples: int) -> np.ndarray:
    edges = np.linspace(0.0, length, n_samples + 1)
    centers = []
    last_center = None
    for i in range(n_samples):
        local = points[(progress >= edges[i]) & (progress <= edges[i + 1])]
        if len(local) >= 4:
            center = robust_local_center(local, last_center, keep_fraction=0.55)
            centers.append(center)
            last_center = center
        elif last_center is not None:
            centers.append(last_center.copy())
    return np.array(centers, dtype=float)


def refine_principal_curve(points: np.ndarray, centers: np.ndarray, passes: int, smoothness: float) -> np.ndarray:
    refined = centers.copy()
    for _ in range(passes):
        assigned = nearest_center_indices(points, refined)
        updated = refined.copy()
        for i in range(1, len(refined) - 1):
            local = points[assigned == i]
            if len(local) >= 4:
                reference = 0.5 * (refined[i - 1] + refined[i + 1])
                data_center = robust_local_center(local, reference, keep_fraction=0.6)
            else:
                data_center = refined[i]
            neighbor_center = 0.5 * (refined[i - 1] + refined[i + 1])
            updated[i] = (1.0 - smoothness) * data_center + smoothness * neighbor_center
        updated[0] = refined[0]
        updated[-1] = refined[-1]
        refined = updated
    return refined


def nearest_center_indices(points: np.ndarray, centers: np.ndarray) -> np.ndarray:
    diff = points[:, None, :] - centers[None, :, :]
    distances = np.einsum("ijk,ijk->ij", diff, diff)
    return np.argmin(distances, axis=1)


def trimmed_mean(points: np.ndarray, trim_fraction: float) -> np.ndarray:
    center = points.mean(axis=0)
    distances = np.linalg.norm(points - center, axis=1)
    keep = max(4, int(len(points) * (1.0 - trim_fraction)))
    indices = np.argsort(distances)[:keep]
    return points[indices].mean(axis=0)


def robust_local_center(points: np.ndarray, reference, keep_fraction: float = 0.6) -> np.ndarray:
    if len(points) == 0:
        raise ValueError("Cannot compute a center from no points.")
    if len(points) < 5:
        return points.mean(axis=0)
    ref = np.array(reference, dtype=float) if reference is not None else np.median(points, axis=0)
    distances = np.linalg.norm(points - ref, axis=1)
    keep = int(np.clip(round(len(points) * keep_fraction), 4, len(points)))
    local = points[np.argsort(distances)[:keep]]
    return trimmed_mean(local, trim_fraction=0.15)


def plane_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = np.array(normal, dtype=float)
    n_norm = np.linalg.norm(n)
    if n_norm == 0:
        n = np.array([1.0, 0.0, 0.0])
    else:
        n /= n_norm
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(n, ref)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    v /= np.linalg.norm(v)
    return u, v


def project_to_plane(points: np.ndarray, normal: np.ndarray) -> np.ndarray:
    origin = points.mean(axis=0)
    u, v = plane_basis(normal)
    shifted = points - origin
    return np.column_stack((shifted @ u, shifted @ v))


def polygon_area_2d(points_2d: np.ndarray) -> float:
    if len(points_2d) < 3:
        return 0.0
    x = points_2d[:, 0]
    y = points_2d[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polygon_centroid_3d(points: np.ndarray, normal: np.ndarray) -> np.ndarray:
    if len(points) < 3:
        return points.mean(axis=0)
    origin = points.mean(axis=0)
    u, v = plane_basis(normal)
    shifted = points - origin
    points_2d = np.column_stack((shifted @ u, shifted @ v))
    area = polygon_area_2d(points_2d)
    if abs(area) < 1e-10:
        return trimmed_mean(points, trim_fraction=0.05)
    x = points_2d[:, 0]
    y = points_2d[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    cx = float(np.sum((x + np.roll(x, -1)) * cross) / (6.0 * area))
    cy = float(np.sum((y + np.roll(y, -1)) * cross) / (6.0 * area))
    return origin + cx * u + cy * v


def generate_surface_path_centerline(surface, source_point, target_point, config: CenterlineConfig) -> vtk.vtkPolyData:
    # Debug-only fallback. This follows the surface and is not a lumen centerline.
    cleaned = clean_surface(surface)
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(cleaned)
    locator.BuildLocator()
    source_id = locator.FindClosestPoint(source_point)
    target_id = locator.FindClosestPoint(target_point)

    path = vtk.vtkDijkstraGraphGeodesicPath()
    path.SetInputData(cleaned)
    path.SetStartVertex(source_id)
    path.SetEndVertex(target_id)
    path.Update()
    polydata = path.GetOutput()
    if polydata is None or polydata.GetNumberOfPoints() == 0:
        raise ValueError("VTK surface-path centerline failed. Check selected points and surface connectivity.")
    return finish_centerline(smooth_centerline(polydata, config.spline_filter_length), config, source_point)


def effective_resampling_step(config: CenterlineConfig, fallback: float = 1e-6) -> float:
    return max(float(config.resampling_step_length), fallback, 1e-6)


def effective_section_radius(surface: vtk.vtkPolyData, config: CenterlineConfig) -> float:
    if config.cross_section_radius > 0:
        return float(config.cross_section_radius)
    return max(surface.GetLength() * 0.08, 1e-6)


def finish_centerline(polydata: vtk.vtkPolyData, config: CenterlineConfig, source_point) -> vtk.vtkPolyData:
    if config.resampling:
        resampled = resample_polyline(polydata, config.resampling_step_length)
        if resampled.GetNumberOfPoints() >= 2:
            polydata = resampled
    return orient_centerline(polydata, source_point)


def enforce_monotone_selected_path(polydata: vtk.vtkPolyData, source_point, target_point, tolerance: float) -> vtk.vtkPolyData:
    points = vtk_points_to_array(polydata.GetPoints())
    if len(points) < 4:
        return polydata
    source = np.array(source_point, dtype=float)
    target = np.array(target_point, dtype=float)
    axis = target - source
    length = float(np.linalg.norm(axis))
    if length <= 0:
        return polydata
    axis /= length
    progress = (points - source) @ axis
    keep = np.ones(len(points), dtype=bool)
    max_progress = progress[0]
    for i in range(1, len(points)):
        if progress[i] + tolerance < max_progress:
            keep[i] = False
        else:
            max_progress = max(max_progress, progress[i])
    filtered = points[keep]
    if len(filtered) < 3 or len(filtered) < 0.6 * len(points):
        return polydata
    return polyline_from_points(filtered)


def vtk_points_to_array(points: vtk.vtkPoints) -> np.ndarray:
    return np.array([points.GetPoint(i) for i in range(points.GetNumberOfPoints())], dtype=float)


def smooth_point_sequence(points: np.ndarray, passes: int = 3, window: int = 5) -> np.ndarray:
    if len(points) <= 4:
        return points
    window = max(3, window | 1)
    radius = window // 2
    smoothed = points.copy()
    for _ in range(passes):
        next_points = smoothed.copy()
        for i in range(1, len(points) - 1):
            lo = max(0, i - radius)
            hi = min(len(points), i + radius + 1)
            next_points[i] = smoothed[lo:hi].mean(axis=0)
        smoothed = next_points
    smoothed[0] = points[0]
    smoothed[-1] = points[-1]
    return smoothed


def polyline_from_points(points: np.ndarray) -> vtk.vtkPolyData:
    vtk_points = vtk.vtkPoints()
    line = vtk.vtkPolyLine()
    line.GetPointIds().SetNumberOfIds(len(points))
    for i, point in enumerate(points):
        vtk_points.InsertNextPoint(float(point[0]), float(point[1]), float(point[2]))
        line.GetPointIds().SetId(i, i)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(line)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetLines(cells)
    return polydata


def polyline_length(polydata: vtk.vtkPolyData) -> float:
    points = polydata.GetPoints()
    if points is None or points.GetNumberOfPoints() < 2:
        return 0.0
    return sum(
        float(np.linalg.norm(np.array(points.GetPoint(i)) - np.array(points.GetPoint(i - 1))))
        for i in range(1, points.GetNumberOfPoints())
    )


def resample_polyline(polydata: vtk.vtkPolyData, step_length: float) -> vtk.vtkPolyData:
    if step_length <= 0:
        return polydata
    spline_filter = vtk.vtkSplineFilter()
    spline_filter.SetInputData(polydata)
    spline_filter.SetSubdivideToLength()
    spline_filter.SetLength(step_length)
    spline_filter.Update()
    return spline_filter.GetOutput()


def smooth_centerline(centerline: vtk.vtkPolyData, length: float) -> vtk.vtkPolyData:
    if length <= 0:
        return centerline
    spline_filter = vtk.vtkSplineFilter()
    spline_filter.SetInputData(centerline)
    spline_filter.SetSubdivideToLength()
    spline_filter.SetLength(length)
    spline_filter.Update()
    return spline_filter.GetOutput()


def orient_centerline(centerline: vtk.vtkPolyData, source_point) -> vtk.vtkPolyData:
    points = centerline.GetPoints()
    n_points = points.GetNumberOfPoints()
    if n_points < 2:
        return centerline

    source = np.array(source_point)
    first = np.array(points.GetPoint(0))
    last = np.array(points.GetPoint(n_points - 1))
    if np.linalg.norm(first - source) <= np.linalg.norm(last - source):
        return centerline

    reversed_points = vtk.vtkPoints()
    reversed_lines = vtk.vtkCellArray()
    for i in range(n_points):
        reversed_points.InsertNextPoint(points.GetPoint(n_points - i - 1))

    line = vtk.vtkPolyLine()
    line.GetPointIds().SetNumberOfIds(n_points)
    for i in range(n_points):
        line.GetPointIds().SetId(i, i)
    reversed_lines.InsertNextCell(line)

    reversed_centerline = vtk.vtkPolyData()
    reversed_centerline.SetPoints(reversed_points)
    reversed_centerline.SetLines(reversed_lines)
    return reversed_centerline


def write_polydata(polydata: vtk.vtkPolyData, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.Write()


def calculate_measurements(surface, centerline, config: CenterlineConfig, on_section=None) -> MeasurementResult:
    points = centerline.GetPoints()
    n_points = points.GetNumberOfPoints() if points else 0
    if n_points < 3:
        raise ValueError("At least three centerline points are required for cross sections.")

    result = MeasurementResult()
    result.points = [np.array(points.GetPoint(i)) for i in range(n_points)]
    result.distances = cumulative_distances(result.points)
    cleaned_surface = clean_surface(surface)

    for index, point in enumerate(result.points):
        plane = vtk.vtkPlane()
        plane.SetOrigin(point)
        plane.SetNormal(tangent_at(result.points, index))
        section = cut_cross_section(cleaned_surface, plane, point, config)
        if section.GetNumberOfPoints() == 0 or section.GetNumberOfPolys() == 0:
            area = 0.0
            perimeter = 0.0
        else:
            area = surface_area(section)
            perimeter = boundary_perimeter(section)
        result.cross_sections.append(section)
        result.cross_section_areas.append(area)
        result.cross_section_perimeters.append(perimeter)
        if on_section is not None:
            on_section(index, section, area, perimeter)

    result.surface_areas = projected_surface_areas(result.distances, result.cross_section_perimeters)
    return result


def cumulative_distances(points: list[np.ndarray]) -> list[float]:
    distances = [0.0]
    for i in range(1, len(points)):
        distances.append(distances[-1] + float(np.linalg.norm(points[i] - points[i - 1])))
    return distances


def tangent_at(points: list[np.ndarray], index: int) -> np.ndarray:
    if index == 0:
        tangent = points[1] - points[0]
    elif index == len(points) - 1:
        tangent = points[-1] - points[-2]
    else:
        tangent = points[index + 1] - points[index - 1]
    norm = np.linalg.norm(tangent)
    return tangent / norm if norm else np.array([1.0, 0.0, 0.0])


def cut_cross_section(surface, plane, point, config: CenterlineConfig) -> vtk.vtkPolyData:
    sphere = vtk.vtkSphere()
    sphere.SetCenter(point)
    sphere.SetRadius(effective_section_radius(surface, config))

    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(surface)
    clipper.SetClipFunction(sphere)
    if config.sphere_inside_out:
        clipper.InsideOutOn()
    else:
        clipper.InsideOutOff()
    clipper.Update()

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(clipper.GetOutput())
    cutter.Update()
    raw_section = cutter.GetOutput()
    if raw_section.GetNumberOfPoints() < 3:
        return vtk.vtkPolyData()

    stripper = vtk.vtkStripper()
    stripper.SetInputData(raw_section)
    stripper.JoinContiguousSegmentsOn()
    stripper.Update()
    if stripper.GetOutput().GetNumberOfPoints() < 3:
        return vtk.vtkPolyData()

    triangulator = vtk.vtkContourTriangulator()
    triangulator.SetInputData(stripper.GetOutput())
    triangulator.Update()
    return triangulator.GetOutput()


def surface_area(polydata: vtk.vtkPolyData) -> float:
    if polydata.GetNumberOfPolys() == 0:
        return 0.0
    mass = vtk.vtkMassProperties()
    mass.SetInputData(polydata)
    return float(mass.GetSurfaceArea())


def boundary_perimeter(polydata: vtk.vtkPolyData) -> float:
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(polydata)
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.BoundaryEdgesOn()
    edges.Update()
    boundary = edges.GetOutput()
    if boundary.GetNumberOfCells() == 0:
        return 0.0

    total = 0.0
    ids = vtk.vtkIdList()
    for i in range(boundary.GetNumberOfCells()):
        boundary.GetCellPoints(i, ids)
        for j in range(ids.GetNumberOfIds() - 1):
            p0 = np.array(boundary.GetPoint(ids.GetId(j)))
            p1 = np.array(boundary.GetPoint(ids.GetId(j + 1)))
            total += float(np.linalg.norm(p1 - p0))
    return total


def projected_surface_areas(distances: list[float], perimeters: list[float]) -> list[float]:
    areas = []
    for i, perimeter in enumerate(perimeters):
        if i == 0:
            areas.append(0.0)
        else:
            distance = distances[i] - distances[i - 1]
            areas.append(0.5 * (perimeter + perimeters[i - 1]) * distance)
    return areas


def save_measurements_csv(result: MeasurementResult, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "distance along curve",
        "x",
        "y",
        "z",
        "cross-sectional area at vertex",
        "cross-section perimeter at vertex",
        "surface area segment projected to this vertex",
    ]
    with open(path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for i, point in enumerate(result.points):
            writer.writerow([
                result.distances[i],
                point[0],
                point[1],
                point[2],
                result.cross_section_areas[i],
                result.cross_section_perimeters[i],
                result.surface_areas[i],
            ])
