from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


Color = tuple[float, float, float]


def normalize_color(value: Iterable[float]) -> Color:
    color = tuple(float(c) for c in value)
    if len(color) != 3:
        raise ValueError("Color must have exactly three channels.")
    if any(c > 1.0 for c in color):
        color = tuple(c / 255.0 for c in color)
    return color  # type: ignore[return-value]


@dataclass(slots=True)
class DisplayConfig:
    background_color: Color = field(default_factory=lambda: (1.0, 1.0, 1.0))
    surface_color: Color = field(default_factory=lambda: normalize_color((157, 195, 231)))
    centerline_color: Color = field(default_factory=lambda: normalize_color((239, 122, 109)))
    selected_point_color: Color = field(default_factory=lambda: normalize_color((215, 99, 160)))
    cross_section_color: Color = field(default_factory=lambda: normalize_color((147, 148, 231)))
    surface_opacity: float = 1.0
    completed_surface_opacity: float = 0.3
    centerline_width: float = 4.0
    cross_section_line_width: float = 2.0
    cross_section_opacity: float = 0.4
    point_radius: float = 0.5


@dataclass(slots=True)
class CenterlineConfig:
    surface_file: Path
    centerline_file: Path
    csv_file: Path | None = None
    centerline_method: str = "section_centroid"
    picker_tolerance: float = 0.005
    resampling_step_length: float = 0.05
    spline_filter_length: float = 0.5
    cross_section_radius: float = 20.0
    seed_selector_name: str = "pointlist"
    append_end_points: int = 1
    resampling: int = 1
    sphere_inside_out: bool = True
    display: DisplayConfig = field(default_factory=DisplayConfig)

    @classmethod
    def from_legacy_output(cls, surface_file: str | Path, output_file: str | Path, **kwargs) -> "CenterlineConfig":
        output = Path(output_file)
        if output.suffix.lower() in {".csv", ".scv"}:
            csv_file = output.with_suffix(".csv")
            centerline_file = output.with_name(f"{output.stem}_centerline.vtk")
        else:
            centerline_file = output
            csv_file = output.with_suffix(".csv")

        display = DisplayConfig(
            background_color=normalize_color(kwargs.pop("interface_bg_color", kwargs.pop("background_color", (255, 255, 255)))),
            surface_color=normalize_color(kwargs.pop("surface_color", (157, 195, 231))),
            centerline_color=normalize_color(kwargs.pop("centerline_color", (239, 122, 109))),
            selected_point_color=normalize_color(kwargs.pop("selected_point_color", (215, 99, 160))),
            cross_section_color=normalize_color(kwargs.pop("cross_section_display_color", (147, 148, 231))),
            cross_section_line_width=float(kwargs.pop("cross_section_line_width", 2.0)),
            cross_section_opacity=float(kwargs.pop("cross_section_opacity", 0.4)),
        )
        return cls(
            surface_file=Path(surface_file),
            centerline_file=centerline_file,
            csv_file=csv_file,
            centerline_method=str(kwargs.pop("centerline_method", "section_centroid")),
            picker_tolerance=float(kwargs.pop("picker_tolerance", 0.005)),
            resampling_step_length=float(kwargs.pop("resampling_step_length", 0.05)),
            spline_filter_length=float(kwargs.pop("spline_filter_length", 0.5)),
            cross_section_radius=float(kwargs.pop("cross_section_radius", 20.0)),
            seed_selector_name=str(kwargs.pop("seed_selector_name", "pointlist")),
            append_end_points=int(kwargs.pop("append_end_points", 1)),
            resampling=int(kwargs.pop("resampling", 1)),
            sphere_inside_out=bool(kwargs.pop("sphere_inside_out", True)),
            display=display,
        )
