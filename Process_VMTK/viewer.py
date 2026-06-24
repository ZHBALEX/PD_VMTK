from __future__ import annotations

import traceback

import vtk

from .config import CenterlineConfig
from .core import (
    MeasurementResult,
    calculate_measurements,
    generate_centerline,
    load_stl,
    save_measurements_csv,
    write_polydata,
)


class CenterlineExtraction:
    """Interactive VTK workflow for point picking, centerlines, and measurements."""

    def __init__(self, surface_file, output_file=None, **kwargs):
        if isinstance(surface_file, CenterlineConfig):
            self.config = surface_file
        else:
            if output_file is None:
                raise ValueError("output_file is required when surface_file is a path.")
            self.config = CenterlineConfig.from_legacy_output(surface_file, output_file, **kwargs)

        self.surface = None
        self.centerline = None
        self.measurements: MeasurementResult | None = None
        self.selected_points = [None, None]
        self.point_actors = [None, None]
        self.cross_section_actors = []
        self.pick_mode = None
        self.pick_observer_id = None

        self.renderer = None
        self.render_window = None
        self.interactor = None
        self.surface_actor = None
        self.centerline_actor = None
        self.picker = vtk.vtkCellPicker()
        self.picker.SetTolerance(self.config.picker_tolerance)
        self.config_provider = None

    def log(self, message: str) -> None:
        print(message)

    def load_surface(self) -> None:
        self.surface = load_stl(self.config.surface_file)
        self.log(f"Loaded surface: {self.config.surface_file}")

    def setup_render(self, render_window=None, interactor=None) -> None:
        if self.surface is None:
            raise RuntimeError("Load a surface before setting up the renderer.")

        display = self.config.display
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.surface)
        self.surface_actor = vtk.vtkActor()
        self.surface_actor.SetMapper(mapper)
        self.surface_actor.GetProperty().SetColor(*display.surface_color)
        self.surface_actor.GetProperty().SetOpacity(display.surface_opacity)

        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(*display.background_color)
        self.renderer.AddActor(self.surface_actor)

        self.render_window = render_window or vtk.vtkRenderWindow()
        self.render_window.GetRenderers().RemoveAllItems()
        self.render_window.AddRenderer(self.renderer)
        self.render_window.SetWindowName("PD VMTK Centerline")
        if render_window is None:
            self.render_window.SetSize(1200, 850)

        self.interactor = interactor or vtk.vtkRenderWindowInteractor()
        if interactor is None:
            self.interactor.SetRenderWindow(self.render_window)
        self.interactor.SetPicker(self.picker)
        self.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        if self.pick_observer_id is not None:
            self.interactor.RemoveObserver(self.pick_observer_id)
        self.pick_observer_id = self.interactor.AddObserver("LeftButtonPressEvent", self.left_button_press_event, 1.0)
        self.renderer.ResetCamera()

    def update_config(self, config: CenterlineConfig) -> None:
        self.config = config
        self.picker.SetTolerance(config.picker_tolerance)
        display = config.display
        if self.renderer is not None:
            self.renderer.SetBackground(*display.background_color)
        if self.surface_actor is not None:
            self.surface_actor.GetProperty().SetColor(*display.surface_color)
        if self.centerline_actor is not None:
            self.centerline_actor.GetProperty().SetColor(*display.centerline_color)
            self.centerline_actor.GetProperty().SetLineWidth(display.centerline_width)
        for actor in self.point_actors:
            if actor is not None:
                actor.GetProperty().SetColor(*display.selected_point_color)
        for actor in self.cross_section_actors:
            actor.GetProperty().SetColor(*display.cross_section_color)
            actor.GetProperty().SetLineWidth(display.cross_section_line_width)
            actor.GetProperty().SetOpacity(display.cross_section_opacity)
        if self.render_window is not None:
            self.render_window.Render()

    def detach(self) -> None:
        if self.interactor is not None and self.pick_observer_id is not None:
            self.interactor.RemoveObserver(self.pick_observer_id)
            self.pick_observer_id = None

    def left_button_press_event(self, obj, event) -> None:
        if self.pick_mode is None:
            return
        click_pos = self.interactor.GetEventPosition()
        self.picker.Pick(click_pos[0], click_pos[1], 0, self.renderer)
        if self.picker.GetCellId() < 0:
            self.log("No valid surface point was selected.")
            return

        point = tuple(float(v) for v in self.picker.GetPickPosition())
        index = 0 if self.pick_mode == "source" else 1
        self.set_selected_point(index, point)
        self.log(f"Selected {self.pick_mode}: {point}")
        self.pick_mode = None

    def sync_config_from_provider(self) -> None:
        if self.config_provider is None:
            return
        try:
            self.update_config(self.config_provider())
        except Exception as exc:
            self.log(f"Could not apply current GUI settings: {exc}")

    def begin_pick(self, mode: str) -> None:
        if mode not in {"source", "target"}:
            raise ValueError(f"Unknown pick mode: {mode}")
        self.sync_config_from_provider()
        self.pick_mode = mode
        self.log(f"Pick mode: click one surface point for {mode}.")

    def set_selected_point(self, index: int, position) -> None:
        old_actor = self.point_actors[index]
        if old_actor is not None:
            self.renderer.RemoveActor(old_actor)
        self.selected_points[index] = position
        self.point_actors[index] = self.create_point_actor(position)
        self.renderer.AddActor(self.point_actors[index])
        self.render_window.Render()

    def create_point_actor(self, position):
        display = self.config.display
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(position)
        sphere.SetRadius(display.point_radius)
        sphere.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*display.selected_point_color)
        return actor

    def undo_point(self) -> None:
        for index in (1, 0):
            actor = self.point_actors[index]
            if actor is not None:
                self.selected_points[index] = None
                self.point_actors[index] = None
                self.renderer.RemoveActor(actor)
                self.render_window.Render()
                self.log("Removed the last selected point.")
                return
        self.log("No selected point to undo.")

    def reset_points(self) -> None:
        self.selected_points = [None, None]
        for actor in self.point_actors:
            if actor is not None:
                self.renderer.RemoveActor(actor)
        self.point_actors = [None, None]
        self.render_window.Render()
        self.log("Selection reset.")

    def generate_centerline(self) -> None:
        try:
            if self.selected_points[0] is None or self.selected_points[1] is None:
                self.log("Select source and target points first.")
                return
            resampling_state = "on" if self.config.resampling else "off"
            section_radius = self.config.cross_section_radius if self.config.cross_section_radius > 0 else "auto"
            self.log(
                "Generating centerline: "
                f"method={self.config.centerline_method}, "
                f"picker_tol={self.config.picker_tolerance:g}, "
                f"resample={resampling_state}, "
                f"step={self.config.resampling_step_length:g}, "
                f"spline={self.config.spline_filter_length:g}, "
                f"section_radius={section_radius}"
            )
            self.clear_cross_sections()
            self.measurements = None
            self.centerline = generate_centerline(self.surface, self.selected_points[0], self.selected_points[1], self.config)
            self.display_centerline()
            self.surface_actor.GetProperty().SetOpacity(self.config.display.completed_surface_opacity)
            self.render_window.Render()
            self.log("Centerline generated. Use Save Centerline to write the VTK file.")
        except Exception as exc:
            self.log(f"Error generating centerline: {exc}")
            traceback.print_exc()

    def save_centerline(self) -> None:
        if self.centerline is None:
            self.log("Generate the centerline before saving it.")
            return
        write_polydata(self.centerline, self.config.centerline_file)
        self.log(f"Centerline saved: {self.config.centerline_file}")

    def display_centerline(self) -> None:
        if self.centerline_actor is not None:
            self.renderer.RemoveActor(self.centerline_actor)

        display = self.config.display
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.centerline)
        self.centerline_actor = vtk.vtkActor()
        self.centerline_actor.SetMapper(mapper)
        self.centerline_actor.GetProperty().SetColor(*display.centerline_color)
        self.centerline_actor.GetProperty().SetLineWidth(display.centerline_width)
        self.renderer.AddActor(self.centerline_actor)

    def calculate_cross_sections(self) -> None:
        try:
            if self.centerline is None:
                self.log("Generate the centerline before calculating cross sections.")
                return
            self.clear_cross_sections()
            self.measurements = calculate_measurements(
                self.surface,
                self.centerline,
                self.config,
                on_section=self._display_section_progress,
            )
            self.log(f"Cross sections calculated: {len(self.measurements.points)}")
        except Exception as exc:
            self.log(f"Error calculating cross sections: {exc}")
            traceback.print_exc()

    def _display_section_progress(self, index, section, area, perimeter) -> None:
        if section.GetNumberOfPoints() > 0 and section.GetNumberOfPolys() > 0:
            self.display_cross_section(section)
        self.log(f"Section {index}: area={area:.6g}, perimeter={perimeter:.6g}")

    def display_cross_section(self, section) -> None:
        display = self.config.display
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(section)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*display.cross_section_color)
        actor.GetProperty().SetLineWidth(display.cross_section_line_width)
        actor.GetProperty().SetOpacity(display.cross_section_opacity)
        actor.GetProperty().SetRepresentationToSurface()
        self.renderer.AddActor(actor)
        self.cross_section_actors.append(actor)
        self.render_window.Render()

    def clear_cross_sections(self) -> None:
        for actor in self.cross_section_actors:
            self.renderer.RemoveActor(actor)
        self.cross_section_actors.clear()

    def save_results(self) -> None:
        if self.measurements is None:
            self.log("Calculate cross sections before saving CSV.")
            return
        if self.config.csv_file is None:
            self.log("No CSV output path configured.")
            return
        save_measurements_csv(self.measurements, self.config.csv_file)
        self.log(f"CSV saved: {self.config.csv_file}")

    def run(self) -> None:
        self.load_surface()
        self.setup_render()
        self.log("Mouse interaction stays in camera mode. Use GUI buttons to pick points, generate, and save.")
        self.log(f"Centerline method: {self.config.centerline_method}")
        self.render_window.Render()
        self.interactor.Start()
