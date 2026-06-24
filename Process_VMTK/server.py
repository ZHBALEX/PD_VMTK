from __future__ import annotations

import json
import mimetypes
import tempfile
import threading
import webbrowser
import argparse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

import vtk

from .config import CenterlineConfig
from .core import (
    calculate_measurements,
    generate_centerline,
    load_stl,
    save_measurements_csv,
    write_polydata,
)


ROOT = Path(__file__).resolve().parent
WEB_ROOT = ROOT / "web"


class AppState:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.surface = None
        self.centerline = None
        self.measurements = None
        self.model_name = "centerline"
        self.surface_path: Path | None = None


STATE = AppState()


class VmtkWebHandler(BaseHTTPRequestHandler):
    server_version = "PDVMTKWeb/1.0"

    def log_message(self, format: str, *args) -> None:
        print(f"{self.address_string()} - {format % args}")

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        path = parsed.path
        if path == "/":
            self._send_file(WEB_ROOT / "index.html")
            return
        if path.startswith("/web/"):
            self._send_file(WEB_ROOT / path[len("/web/"):])
            return
        if path == "/api/health":
            self._send_json({"ok": True})
            return
        if path == "/api/surface.stl":
            self._send_current_surface()
            return
        self._send_json({"error": "Not found"}, status=404)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        try:
            if parsed.path == "/api/model":
                self._load_model()
            elif parsed.path == "/api/centerline":
                self._centerline()
            elif parsed.path == "/api/save-centerline":
                self._save_centerline()
            elif parsed.path == "/api/cross-sections":
                self._cross_sections()
            elif parsed.path == "/api/save-csv":
                self._save_csv()
            elif parsed.path == "/api/cut-box":
                self._cut_box()
            elif parsed.path == "/api/transform":
                self._transform_surface()
            elif parsed.path == "/api/save-stl":
                self._save_stl()
            else:
                self._send_json({"error": "Not found"}, status=404)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=500)

    def _load_model(self) -> None:
        length = int(self.headers.get("Content-Length", "0"))
        if length <= 0:
            raise ValueError("No STL data was received.")
        raw_name = unquote(self.headers.get("X-Filename", "model.stl"))
        model_name = Path(raw_name).stem or "centerline"
        data = self.rfile.read(length)

        temp = tempfile.NamedTemporaryFile(delete=False, suffix=".stl", prefix="pd_vmtk_")
        try:
            temp.write(data)
            temp_path = Path(temp.name)
        finally:
            temp.close()

        surface = load_stl(temp_path)
        with STATE.lock:
            STATE.surface = surface
            STATE.centerline = None
            STATE.measurements = None
            STATE.model_name = model_name
            STATE.surface_path = temp_path

        self._send_json({
            "ok": True,
            "modelName": model_name,
            "points": surface.GetNumberOfPoints(),
            "cells": surface.GetNumberOfCells(),
            "bounds": list(surface.GetBounds()),
        })

    def _centerline(self) -> None:
        payload = self._read_json()
        source = payload.get("source")
        target = payload.get("target")
        if not source or not target:
            raise ValueError("Pick source and target points before generating centerline.")

        with STATE.lock:
            surface = STATE.surface
            model_name = STATE.model_name
        if surface is None:
            raise ValueError("Load an STL model first.")

        config = self._config_from_payload(payload, model_name)
        centerline = generate_centerline(surface, source, target, config)
        with STATE.lock:
            STATE.centerline = centerline
            STATE.measurements = None

        self._send_json({
            "ok": True,
            "points": polydata_points(centerline),
            "count": centerline.GetNumberOfPoints(),
            "settings": summarize_config(config),
        })

    def _save_centerline(self) -> None:
        payload = self._read_json()
        output_dir = Path(payload.get("outputDir") or Path.cwd())
        output_name = (payload.get("outputName") or STATE.model_name or "centerline").strip()
        if not output_name:
            output_name = "centerline"
        path = output_dir / f"{output_name}_centerline.vtk"

        with STATE.lock:
            centerline = STATE.centerline
        if centerline is None:
            raise ValueError("Generate a centerline before saving.")

        write_polydata(centerline, path)
        self._send_json({"ok": True, "path": str(path)})

    def _cross_sections(self) -> None:
        payload = self._read_json()
        with STATE.lock:
            surface = STATE.surface
            centerline = STATE.centerline
            model_name = STATE.model_name
        if surface is None:
            raise ValueError("Load an STL model first.")
        if centerline is None:
            raise ValueError("Generate a centerline before calculating cross sections.")

        config = self._config_from_payload(payload, model_name)
        measurements = calculate_measurements(surface, centerline, config)
        with STATE.lock:
            STATE.measurements = measurements

        self._send_json({
            "ok": True,
            "count": len(measurements.points),
            "areas": measurements.cross_section_areas,
            "perimeters": measurements.cross_section_perimeters,
            "distances": measurements.distances,
            "sections": [section_boundary_polylines(section) for section in measurements.cross_sections],
        })

    def _save_csv(self) -> None:
        payload = self._read_json()
        output_dir = Path(payload.get("outputDir") or Path.cwd())
        output_name = (payload.get("outputName") or STATE.model_name or "centerline").strip() or "centerline"
        path = output_dir / f"{output_name}.csv"

        with STATE.lock:
            measurements = STATE.measurements
        if measurements is None:
            raise ValueError("Calculate cross sections before saving CSV.")

        save_measurements_csv(measurements, path)
        self._send_json({"ok": True, "path": str(path)})

    def _cut_box(self) -> None:
        payload = self._read_json()
        bounds = payload.get("bounds")
        mode = payload.get("mode", "outside")
        resample = bool(payload.get("resample", True))
        if not bounds or len(bounds) != 6:
            raise ValueError("Cut box bounds must contain [xmin, xmax, ymin, ymax, zmin, zmax].")
        bounds = [float(v) for v in bounds]
        if bounds[0] >= bounds[1] or bounds[2] >= bounds[3] or bounds[4] >= bounds[5]:
            raise ValueError("Cut box min values must be smaller than max values.")

        with STATE.lock:
            surface = STATE.surface
        if surface is None:
            raise ValueError("Load an STL model first.")

        cut_surface = cut_surface_with_box(surface, bounds, keep_inside=(mode == "inside"), resample=resample)
        if cut_surface.GetNumberOfPoints() == 0:
            raise ValueError("Cut produced an empty surface. Check the box range or cut mode.")

        with STATE.lock:
            STATE.surface = cut_surface
            STATE.centerline = None
            STATE.measurements = None

        self._send_json({
            "ok": True,
            "points": cut_surface.GetNumberOfPoints(),
            "cells": cut_surface.GetNumberOfCells(),
            "bounds": list(cut_surface.GetBounds()),
        })

    def _transform_surface(self) -> None:
        payload = self._read_json()
        with STATE.lock:
            surface = STATE.surface
        if surface is None:
            raise ValueError("Load an STL model first.")

        action = payload.get("action", "apply")
        if action == "origin":
            bounds = surface.GetBounds()
            center = (
                0.5 * (bounds[0] + bounds[1]),
                0.5 * (bounds[2] + bounds[3]),
                0.5 * (bounds[4] + bounds[5]),
            )
            transform = vtk.vtkTransform()
            transform.Translate(-center[0], -center[1], -center[2])
        elif action == "center":
            target = [float(v) for v in payload.get("translate", [0, 0, 0])]
            bounds = surface.GetBounds()
            center = (
                0.5 * (bounds[0] + bounds[1]),
                0.5 * (bounds[2] + bounds[3]),
                0.5 * (bounds[4] + bounds[5]),
            )
            transform = vtk.vtkTransform()
            transform.Translate(target[0] - center[0], target[1] - center[1], target[2] - center[2])
        else:
            translate = [float(v) for v in payload.get("translate", [0, 0, 0])]
            rotate = [float(v) for v in payload.get("rotate", [0, 0, 0])]
            scale = [float(v) for v in payload.get("scale", [1, 1, 1])]
            bounds = surface.GetBounds()
            center = (
                0.5 * (bounds[0] + bounds[1]),
                0.5 * (bounds[2] + bounds[3]),
                0.5 * (bounds[4] + bounds[5]),
            )
            transform = vtk.vtkTransform()
            transform.PostMultiply()
            transform.Translate(-center[0], -center[1], -center[2])
            transform.Scale(scale)
            transform.RotateX(rotate[0])
            transform.RotateY(rotate[1])
            transform.RotateZ(rotate[2])
            transform.Translate(center[0] + translate[0], center[1] + translate[1], center[2] + translate[2])

        transformed = transform_polydata(surface, transform)
        with STATE.lock:
            STATE.surface = transformed
            STATE.centerline = None
            STATE.measurements = None
        self._send_json({
            "ok": True,
            "points": transformed.GetNumberOfPoints(),
            "cells": transformed.GetNumberOfCells(),
            "bounds": list(transformed.GetBounds()),
        })

    def _save_stl(self) -> None:
        payload = self._read_json()
        output_dir = Path(payload.get("outputDir") or Path.cwd())
        output_name = (payload.get("outputName") or STATE.model_name or "geometry").strip() or "geometry"
        path = output_dir / f"{output_name}_edited.stl"
        with STATE.lock:
            surface = STATE.surface
        if surface is None:
            raise ValueError("Load an STL model before saving geometry.")
        write_stl(surface, path)
        self._send_json({"ok": True, "path": str(path)})

    def _send_current_surface(self) -> None:
        with STATE.lock:
            surface = STATE.surface
        if surface is None:
            self._send_json({"error": "No surface loaded."}, status=404)
            return
        data = stl_bytes(surface)
        self.send_response(200)
        self.send_header("Content-Type", "model/stl")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _config_from_payload(self, payload: dict, model_name: str) -> CenterlineConfig:
        config = payload.get("config") or {}
        output_dir = Path(payload.get("outputDir") or Path.cwd())
        output_name = (payload.get("outputName") or model_name or "centerline").strip() or "centerline"
        return CenterlineConfig(
            surface_file=Path("uploaded.stl"),
            centerline_file=output_dir / f"{output_name}_centerline.vtk",
            csv_file=output_dir / f"{output_name}.csv",
            centerline_method=str(config.get("method", "section_centroid")),
            picker_tolerance=float(config.get("pickerTolerance", 0.005)),
            resampling_step_length=float(config.get("resamplingStep", 0.05)),
            spline_filter_length=float(config.get("splineLength", 0.5)),
            cross_section_radius=float(config.get("sectionRadius", 20.0)),
            append_end_points=1 if config.get("appendEndPoints", True) else 0,
            resampling=1 if config.get("resample", True) else 0,
            sphere_inside_out=bool(config.get("clipInside", True)),
        )

    def _read_json(self) -> dict:
        length = int(self.headers.get("Content-Length", "0"))
        if length <= 0:
            return {}
        return json.loads(self.rfile.read(length).decode("utf-8"))

    def _send_file(self, path: Path) -> None:
        resolved = path.resolve()
        if not str(resolved).startswith(str(WEB_ROOT.resolve())) or not resolved.exists() or not resolved.is_file():
            self._send_json({"error": "Not found"}, status=404)
            return
        content_type = mimetypes.guess_type(resolved.name)[0] or "application/octet-stream"
        data = resolved.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _send_json(self, payload: dict, status: int = 200) -> None:
        data = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)


def polydata_points(polydata: vtk.vtkPolyData) -> list[list[float]]:
    points = polydata.GetPoints()
    if points is None:
        return []
    return [[float(v) for v in points.GetPoint(i)] for i in range(points.GetNumberOfPoints())]


def section_boundary_polylines(polydata: vtk.vtkPolyData) -> list[list[list[float]]]:
    if polydata.GetNumberOfPoints() == 0:
        return []
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(polydata)
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.BoundaryEdgesOn()
    edges.Update()
    stripper = vtk.vtkStripper()
    stripper.SetInputData(edges.GetOutput())
    stripper.JoinContiguousSegmentsOn()
    stripper.Update()
    boundary = stripper.GetOutput()
    lines = []
    ids = vtk.vtkIdList()
    for cell_id in range(boundary.GetNumberOfCells()):
        boundary.GetCellPoints(cell_id, ids)
        if ids.GetNumberOfIds() < 2:
            continue
        line = [[float(v) for v in boundary.GetPoint(ids.GetId(i))] for i in range(ids.GetNumberOfIds())]
        lines.append(line)
    return lines


def stl_bytes(surface: vtk.vtkPolyData) -> bytes:
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=".stl", prefix="pd_vmtk_surface_")
    temp_path = Path(temp.name)
    temp.close()
    try:
        writer = vtk.vtkSTLWriter()
        writer.SetFileName(str(temp_path))
        writer.SetInputData(surface)
        writer.SetFileTypeToBinary()
        writer.Write()
        return temp_path.read_bytes()
    finally:
        try:
            temp_path.unlink()
        except OSError:
            pass


def write_stl(surface: vtk.vtkPolyData, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(surface)
    writer.SetFileTypeToBinary()
    writer.Write()


def transform_polydata(surface: vtk.vtkPolyData, transform: vtk.vtkTransform) -> vtk.vtkPolyData:
    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputData(surface)
    transform_filter.SetTransform(transform)
    transform_filter.Update()
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(transform_filter.GetOutput())
    clean.Update()
    return clean.GetOutput()


def cut_surface_with_box(surface: vtk.vtkPolyData, bounds: list[float], keep_inside: bool, resample: bool) -> vtk.vtkPolyData:
    if keep_inside:
        clipped = clip_closed_surface_to_box(surface, bounds)
    else:
        box = vtk.vtkBox()
        box.SetBounds(bounds)
        clipper = vtk.vtkClipPolyData()
        clipper.SetInputData(surface)
        clipper.SetClipFunction(box)
        clipper.SetValue(0.0)
        clipper.InsideOutOff()
        clipper.Update()
        clipped = clipper.GetOutput()

    triangle = vtk.vtkTriangleFilter()
    triangle.SetInputData(clipped)
    triangle.Update()
    polygons = polygon_only_polydata(triangle.GetOutput())

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(polygons)
    clean.Update()
    result = polygon_only_polydata(clean.GetOutput())

    if resample and result.GetNumberOfPoints() < 250000:
        subdivide = vtk.vtkLinearSubdivisionFilter()
        subdivide.SetInputData(result)
        subdivide.SetNumberOfSubdivisions(1)
        subdivide.Update()
        result = subdivide.GetOutput()

    final_clean = vtk.vtkCleanPolyData()
    final_clean.SetInputData(result)
    final_clean.Update()
    return polygon_only_polydata(final_clean.GetOutput())


def polygon_only_polydata(polydata: vtk.vtkPolyData) -> vtk.vtkPolyData:
    output = vtk.vtkPolyData()
    output.SetPoints(polydata.GetPoints())
    output.SetPolys(polydata.GetPolys())
    return output


def clip_closed_surface_to_box(surface: vtk.vtkPolyData, bounds: list[float]) -> vtk.vtkPolyData:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    span = max(xmax - xmin, ymax - ymin, zmax - zmin, 1.0)
    eps = span * 1e-7
    plane_specs = [
        ((xmin - eps, 0.0, 0.0), (1.0, 0.0, 0.0)),
        ((xmax + eps, 0.0, 0.0), (-1.0, 0.0, 0.0)),
        ((0.0, ymin - eps, 0.0), (0.0, 1.0, 0.0)),
        ((0.0, ymax + eps, 0.0), (0.0, -1.0, 0.0)),
        ((0.0, 0.0, zmin - eps), (0.0, 0.0, 1.0)),
        ((0.0, 0.0, zmax + eps), (0.0, 0.0, -1.0)),
    ]
    planes = vtk.vtkPlaneCollection()
    for origin, normal in plane_specs:
        plane = vtk.vtkPlane()
        plane.SetOrigin(origin)
        plane.SetNormal(normal)
        planes.AddItem(plane)

    clipper = vtk.vtkClipClosedSurface()
    clipper.SetInputData(surface)
    clipper.SetClippingPlanes(planes)
    clipper.GenerateFacesOn()
    clipper.Update()
    return clipper.GetOutput()


def summarize_config(config: CenterlineConfig) -> dict:
    return {
        "method": config.centerline_method,
        "pickerTolerance": config.picker_tolerance,
        "resample": bool(config.resampling),
        "resamplingStep": config.resampling_step_length,
        "splineLength": config.spline_filter_length,
        "sectionRadius": config.cross_section_radius,
    }


def run(host: str = "127.0.0.1", port: int = 8765, open_browser: bool = True) -> None:
    try:
        server = ThreadingHTTPServer((host, port), VmtkWebHandler)
    except OSError:
        server = ThreadingHTTPServer((host, 0), VmtkWebHandler)
    url = f"http://{host}:{server.server_port}/"
    print(f"PD VMTK web app running at {url}")
    if open_browser:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.")
    finally:
        server.server_close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the PD VMTK Three.js web app.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    run(host=args.host, port=args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    main()
