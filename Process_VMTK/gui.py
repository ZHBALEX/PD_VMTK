from __future__ import annotations

import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import vtk
try:
    from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
except Exception:
    vtkTkRenderWindowInteractor = None

from .config import CenterlineConfig, DisplayConfig, normalize_color
from .viewer import CenterlineExtraction


class NativeVtkFrame(tk.Frame):
    """Fallback VTK host for wheels without vtkRenderingTk.dll."""

    def __init__(self, master, width=900, height=620):
        super().__init__(master, width=width, height=height, bg="#eef1ef")
        self.render_window = vtk.vtkRenderWindow()
        self.interactor = vtk.vtkRenderWindowInteractor()
        self._initialized = False
        self.bind("<Map>", self._initialize)
        self.bind("<Configure>", self._resize)

    def _initialize(self, event=None):
        if self._initialized:
            return
        self.update_idletasks()
        hwnd = str(self.winfo_id())
        self.render_window.SetParentInfo(hwnd)
        self.render_window.SetSize(max(self.winfo_width(), 1), max(self.winfo_height(), 1))
        self.interactor.SetRenderWindow(self.render_window)
        self.interactor.Initialize()
        self._initialized = True

    def _resize(self, event):
        if self._initialized:
            self.render_window.SetSize(max(event.width, 1), max(event.height, 1))
            self.render_window.Render()

    def GetRenderWindow(self):
        self._initialize()
        return self.render_window

    def Initialize(self):
        self._initialize()


class VmtkApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PD VMTK Centerline")
        self.geometry("980x680")
        self.minsize(860, 600)
        self.configure(bg="#f7f7f4")
        self._build_state()
        self._build_ui()

    def _build_state(self) -> None:
        self.surface_file = tk.StringVar()
        self.output_dir = tk.StringVar(value=str(Path.cwd()))
        self.output_name = tk.StringVar(value="centerline")
        self.centerline_method = tk.StringVar(value="section_centroid")
        self.picker_tolerance = tk.DoubleVar(value=0.005)
        self.resampling_step_length = tk.DoubleVar(value=0.05)
        self.spline_filter_length = tk.DoubleVar(value=0.5)
        self.cross_section_radius = tk.DoubleVar(value=20.0)
        self.append_end_points = tk.IntVar(value=1)
        self.resampling = tk.IntVar(value=1)
        self.sphere_inside_out = tk.BooleanVar(value=True)
        self.background_color = tk.StringVar(value="255,255,255")
        self.surface_color = tk.StringVar(value="157,195,231")
        self.centerline_color = tk.StringVar(value="239,122,109")
        self.point_color = tk.StringVar(value="215,99,160")
        self.section_color = tk.StringVar(value="147,148,231")
        self.section_opacity = tk.DoubleVar(value=0.4)
        self.extraction = None
        self.vtk_widget = None
        self._startup_messages = []

    def _build_ui(self) -> None:
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("TFrame", background="#f7f7f4")
        style.configure("Side.TFrame", background="#fbfbf8")
        style.configure("Panel.TLabelframe", background="#ffffff", bordercolor="#d8ddd8")
        style.configure("Panel.TLabelframe.Label", background="#ffffff", foreground="#67727e")
        style.configure("TLabel", background="#ffffff", foreground="#67727e")
        style.configure("TButton", padding=7)

        shell = ttk.Frame(self)
        shell.pack(fill="both", expand=True)
        sidebar_canvas = tk.Canvas(shell, width=380, bg="#fbfbf8", highlightthickness=0)
        sidebar_scroll = ttk.Scrollbar(shell, orient="vertical", command=sidebar_canvas.yview)
        sidebar_canvas.configure(yscrollcommand=sidebar_scroll.set)
        sidebar_canvas.pack(side="left", fill="y")
        sidebar_scroll.pack(side="left", fill="y")
        sidebar = ttk.Frame(sidebar_canvas, style="Side.TFrame", padding=18)
        sidebar_window = sidebar_canvas.create_window((0, 0), window=sidebar, anchor="nw")
        sidebar.bind(
            "<Configure>",
            lambda event: sidebar_canvas.configure(scrollregion=sidebar_canvas.bbox("all")),
        )
        sidebar_canvas.bind(
            "<Configure>",
            lambda event: sidebar_canvas.itemconfigure(sidebar_window, width=event.width),
        )
        sidebar_canvas.bind_all("<MouseWheel>", lambda event: sidebar_canvas.yview_scroll(int(-event.delta / 120), "units"))
        workspace = ttk.Frame(shell, padding=18)
        workspace.pack(side="left", fill="both", expand=True)

        header = ttk.Frame(sidebar, style="Side.TFrame")
        header.pack(fill="x", pady=(0, 14))
        tk.Label(header, text="VMTK Centerline", bg="#fbfbf8", fg="#1f2933", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        tk.Label(
            header,
            text="Pick two STL surface points, extract centerline, calculate sections, export CSV.",
            bg="#fbfbf8",
            fg="#67727e",
            justify="left",
            wraplength=330,
        ).pack(anchor="w", pady=(4, 0))

        self._file_panel(sidebar)
        self._parameter_panel(sidebar)
        self._display_panel(sidebar)
        self._action_panel(sidebar)
        self._workspace(workspace)

    def _file_panel(self, parent) -> None:
        panel = self._panel(parent, "Files")
        self._path_row(panel, "STL surface", self.surface_file, self._choose_surface)
        self._path_row(panel, "Output folder", self.output_dir, self._choose_output_dir)
        self._entry(panel, "Output name", self.output_name)

    def _parameter_panel(self, parent) -> None:
        panel = self._panel(parent, "Centerline")
        ttk.Label(panel, text="Method").pack(anchor="w")
        method = ttk.Combobox(
            panel,
            textvariable=self.centerline_method,
            values=("section_centroid", "centroid", "vmtk", "vtk_surface"),
            state="readonly",
        )
        method.pack(fill="x", pady=(0, 8))
        grid = ttk.Frame(panel)
        grid.pack(fill="x")
        self._entry(grid, "Picker tolerance", self.picker_tolerance, 0, 0)
        self._entry(grid, "Resampling step", self.resampling_step_length, 0, 1)
        self._entry(grid, "Spline length", self.spline_filter_length, 1, 0)
        self._entry(grid, "Section radius (0 auto)", self.cross_section_radius, 1, 1)
        ttk.Checkbutton(panel, text="Append end points", variable=self.append_end_points).pack(anchor="w", pady=(6, 0))
        ttk.Checkbutton(panel, text="Resample centerline", variable=self.resampling).pack(anchor="w")
        ttk.Checkbutton(panel, text="Clip inside section sphere", variable=self.sphere_inside_out).pack(anchor="w")

    def _display_panel(self, parent) -> None:
        panel = self._panel(parent, "Display")
        grid = ttk.Frame(panel)
        grid.pack(fill="x")
        self._entry(grid, "Background RGB", self.background_color, 0, 0)
        self._entry(grid, "Surface RGB", self.surface_color, 0, 1)
        self._entry(grid, "Centerline RGB", self.centerline_color, 1, 0)
        self._entry(grid, "Point RGB", self.point_color, 1, 1)
        self._entry(grid, "Section RGB", self.section_color, 2, 0)
        self._entry(grid, "Section opacity", self.section_opacity, 2, 1)

    def _action_panel(self, parent) -> None:
        panel = self._panel(parent, "Run")
        ttk.Button(panel, text="Load / Reload STL", command=self._load_surface_into_viewer).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Pick Source", command=lambda: self._begin_pick("source")).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Pick Target", command=lambda: self._begin_pick("target")).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Generate Centerline", command=self._generate_centerline).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Save Centerline", command=self._save_centerline).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Cross Sections", command=self._calculate_cross_sections).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Save CSV", command=self._save_results).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Reset Points", command=self._reset_points).pack(fill="x", pady=(0, 8))
        ttk.Button(panel, text="Clear Log", command=lambda: self.log_text.delete("1.0", "end")).pack(fill="x")

    def _workspace(self, parent) -> None:
        toolbar = ttk.Frame(parent)
        toolbar.pack(fill="x", pady=(0, 10))
        tk.Label(toolbar, text="3D Viewer", bg="#f7f7f4", fg="#1f2933", font=("Segoe UI", 14, "bold")).pack(side="left")
        tk.Label(toolbar, text="Rotate normally. Use Pick Source / Pick Target buttons for point selection.", bg="#f7f7f4", fg="#67727e").pack(side="left", padx=(12, 0))

        self.viewer_frame = ttk.Frame(parent)
        self.viewer_frame.pack(fill="both", expand=True)
        self.vtk_widget = self._create_vtk_widget(self.viewer_frame)
        self.vtk_widget.pack(fill="both", expand=True)
        self.vtk_widget.Initialize()

        self.log_text = tk.Text(parent, height=7, wrap="word", bg="#17211d", fg="#d8f3dc", insertbackground="#d8f3dc", relief="flat")
        self.log_text.pack(fill="x", pady=(12, 0))
        self._log("Centerline method defaults to section_centroid: local cross-section centers guided by the picked path.")
        for message in self._startup_messages:
            self._log(message)
        self._show_empty_viewer()

    def _create_vtk_widget(self, parent):
        if vtkTkRenderWindowInteractor is not None:
            try:
                return vtkTkRenderWindowInteractor(parent, rw=vtk.vtkRenderWindow(), width=900, height=620)
            except tk.TclError as exc:
                self._startup_messages.append(f"VTK Tk widget unavailable, using native embed fallback: {exc}")
        return NativeVtkFrame(parent, width=900, height=620)

    def _panel(self, parent, title: str):
        panel = ttk.LabelFrame(parent, text=title, style="Panel.TLabelframe", padding=12)
        panel.pack(fill="x", pady=(0, 12))
        return panel

    def _entry(self, parent, label: str, variable, row=None, column=None):
        frame = ttk.Frame(parent)
        if row is None:
            frame.pack(fill="x", pady=(0, 8))
        else:
            frame.grid(row=row, column=column, sticky="ew", padx=(0 if column == 0 else 8, 0), pady=(0, 8))
            parent.columnconfigure(column, weight=1)
        ttk.Label(frame, text=label).pack(anchor="w")
        ttk.Entry(frame, textvariable=variable).pack(fill="x")

    def _path_row(self, parent, label: str, variable, command) -> None:
        frame = ttk.Frame(parent)
        frame.pack(fill="x", pady=(0, 8))
        ttk.Label(frame, text=label).pack(anchor="w")
        row = ttk.Frame(frame)
        row.pack(fill="x")
        ttk.Entry(row, textvariable=variable).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Browse", command=command, width=9).pack(side="left", padx=(8, 0))

    def _choose_surface(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("STL files", "*.stl"), ("All files", "*.*")])
        if path:
            self.surface_file.set(path)
            self.output_name.set(Path(path).stem)
            self._log(f"Selected STL: {path}")
            self._load_surface_into_viewer()

    def _choose_output_dir(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.output_dir.set(path)

    def _load_surface_into_viewer(self) -> None:
        try:
            self.update_idletasks()
            config = self._read_config()
        except Exception as exc:
            messagebox.showerror("Invalid settings", str(exc))
            return
        try:
            if self.extraction is not None:
                self.extraction.detach()
            self.extraction = CenterlineExtraction(config)
            self.extraction.log = self._log
            self.extraction.config_provider = self._read_config
            self.extraction.load_surface()
            self.extraction.setup_render(
                render_window=self.vtk_widget.GetRenderWindow(),
                interactor=self._vtk_interactor(),
            )
            self.vtk_widget.GetRenderWindow().Render()
            try:
                self.vtk_widget.focus_set()
            except Exception:
                pass
            self._log("3D model loaded in the main window.")
            self._log(self._settings_summary(config))
            self._log("Use mouse to rotate/zoom/pan. Use Pick Source / Pick Target buttons to select endpoints.")
        except Exception as exc:
            self._log(f"Viewer failed: {exc}")

    def _vtk_interactor(self):
        return getattr(self.vtk_widget, "interactor", self.vtk_widget)

    def _show_empty_viewer(self) -> None:
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(0.93, 0.95, 0.94)
        self.vtk_widget.GetRenderWindow().GetRenderers().RemoveAllItems()
        self.vtk_widget.GetRenderWindow().AddRenderer(renderer)
        self.vtk_widget.GetRenderWindow().Render()

    def _begin_pick(self, mode: str) -> None:
        if self._need_model():
            return
        if not self._sync_settings():
            return
        self.extraction.begin_pick(mode)

    def _generate_centerline(self) -> None:
        if self._need_model():
            return
        if not self._sync_settings():
            return
        self.extraction.generate_centerline()

    def _save_centerline(self) -> None:
        if self._need_model():
            return
        if not self._sync_settings():
            return
        self.extraction.save_centerline()

    def _calculate_cross_sections(self) -> None:
        if self._need_model():
            return
        if not self._sync_settings():
            return
        self.extraction.calculate_cross_sections()

    def _save_results(self) -> None:
        if self._need_model():
            return
        if not self._sync_settings():
            return
        self.extraction.save_results()

    def _reset_points(self) -> None:
        if self._need_model():
            return
        self._sync_settings(show_error=False)
        self.extraction.reset_points()

    def _need_model(self) -> bool:
        if self.extraction is None:
            messagebox.showinfo("No model", "Please choose an STL file first.")
            return True
        return False

    def _sync_settings(self, show_error: bool = True) -> bool:
        try:
            self.update_idletasks()
            config = self._read_config()
        except Exception as exc:
            if show_error:
                messagebox.showerror("Invalid settings", str(exc))
            return False
        self.extraction.update_config(config)
        self._log(self._settings_summary(config))
        return True

    def _settings_summary(self, config: CenterlineConfig) -> str:
        resampling_state = "on" if config.resampling else "off"
        radius = f"{config.cross_section_radius:g}" if config.cross_section_radius > 0 else "auto"
        return (
            "Settings applied: "
            f"method={config.centerline_method}, "
            f"picker_tol={config.picker_tolerance:g}, "
            f"resample={resampling_state}, "
            f"step={config.resampling_step_length:g}, "
            f"spline={config.spline_filter_length:g}, "
            f"section_radius={radius}"
        )

    def _read_config(self) -> CenterlineConfig:
        surface = Path(self.surface_file.get())
        if not surface.exists():
            raise ValueError("Please choose an existing STL surface file.")
        output_dir = Path(self.output_dir.get())
        name = self.output_name.get().strip() or surface.stem
        picker_tolerance = float(self.picker_tolerance.get())
        resampling_step = float(self.resampling_step_length.get())
        spline_length = float(self.spline_filter_length.get())
        section_radius = float(self.cross_section_radius.get())
        section_opacity = float(self.section_opacity.get())
        if picker_tolerance <= 0:
            raise ValueError("Picker tolerance must be greater than 0.")
        if resampling_step <= 0:
            raise ValueError("Resampling step must be greater than 0.")
        if spline_length < 0:
            raise ValueError("Spline length must be 0 or greater.")
        if section_radius < 0:
            raise ValueError("Section radius must be 0 or greater.")
        if not 0 <= section_opacity <= 1:
            raise ValueError("Section opacity must be between 0 and 1.")
        display = DisplayConfig(
            background_color=self._rgb(self.background_color.get()),
            surface_color=self._rgb(self.surface_color.get()),
            centerline_color=self._rgb(self.centerline_color.get()),
            selected_point_color=self._rgb(self.point_color.get()),
            cross_section_color=self._rgb(self.section_color.get()),
            cross_section_opacity=section_opacity,
        )
        return CenterlineConfig(
            surface_file=surface,
            centerline_file=output_dir / f"{name}_centerline.vtk",
            csv_file=output_dir / f"{name}.csv",
            centerline_method=self.centerline_method.get(),
            picker_tolerance=picker_tolerance,
            resampling_step_length=resampling_step,
            spline_filter_length=spline_length,
            cross_section_radius=section_radius,
            append_end_points=int(self.append_end_points.get()),
            resampling=int(self.resampling.get()),
            sphere_inside_out=bool(self.sphere_inside_out.get()),
            display=display,
        )

    def _rgb(self, text: str):
        return normalize_color(float(part.strip()) for part in text.split(","))

    def _log(self, message: str) -> None:
        self.log_text.insert("end", message + "\n")
        self.log_text.see("end")


def main() -> None:
    VmtkApp().mainloop()


if __name__ == "__main__":
    main()
