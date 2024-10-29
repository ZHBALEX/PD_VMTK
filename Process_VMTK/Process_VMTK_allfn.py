import vtk
from vmtk import vmtkscripts
import csv
import os
import sys
import traceback
import numpy as np

class CenterlineExtraction:
    def __init__(self, surface_file, output_file):
        self.surface_file = surface_file
        self.output_file = output_file
        self.surface = None
        self.renderer = None
        self.render_window = None
        self.interactor = None
        self.selected_points = []
        self.point_actors = []
        self.picker = vtk.vtkCellPicker()
        self.picker.SetTolerance(0.005)  # 增大容差值
        self.centerline_actor = None
        self.surface_actor = None
        self.centerlines = None  # Store the centerlines object
        self.cross_section_actors = []
        self.cross_section_areas = []
        self.distance_along_curve = []
        self.cross_section_points = []
        self.surface_areas = []

    def load_surface(self):
        try:
            # Read the STL file
            reader = vtk.vtkSTLReader()
            reader.SetFileName(self.surface_file)
            reader.Update()
            self.surface = reader.GetOutput()
            if self.surface.GetNumberOfPoints() == 0:
                raise ValueError("表面模型未能正确加载，请检查文件路径和文件格式。")
        except Exception as e:
            print("加载表面模型时出错：{}".format(e))
            traceback.print_exc()
            # 不退出程序，允许用户继续操作

    def preprocess_surface(self):
        # 使用 vtkCleanPolyData 来清理表面，移除非流形的几何体
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(self.surface)
        cleaner.Update()
        self.surface = cleaner.GetOutput()

    def setup_render(self):
        # Create a mapper and actor for the surface
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.surface)
        self.surface_actor = vtk.vtkActor()
        self.surface_actor.SetMapper(mapper)

        # Create renderer, render window, and interactor
        self.renderer = vtk.vtkRenderer()
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.render_window)

        # Add the surface actor to the scene
        self.renderer.AddActor(self.surface_actor)
        self.renderer.SetBackground(0.1, 0.2, 0.4)

        # Set up the picker
        self.interactor.SetPicker(self.picker)

        # Use trackball camera style for rotation
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        # Set up event handling
        self.interactor.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.interactor.AddObserver("KeyPressEvent", self.key_press_event)

    def left_button_press_event(self, obj, event):
        click_pos = self.interactor.GetEventPosition()
        self.picker.Pick(click_pos[0], click_pos[1], 0, self.renderer)
        world_position = self.picker.GetPickPosition()
        cell_id = self.picker.GetCellId()
        if cell_id >= 0:
            # Add the picked point to the list
            self.selected_points.append(world_position)
            print("选择的点坐标: {}".format(world_position))
            # Display the picked point
            self.display_point(world_position)
        else:
            print("未选取到有效的表面点。")

    def key_press_event(self, obj, event):
        key = self.interactor.GetKeySym()
        if key == 'space':
            if len(self.selected_points) == 2:
                print("已选择两个点，正在生成中心线...")
                self.generate_centerline()
            else:
                print("已选择{}个点。请选择两个点以继续。".format(len(self.selected_points)))
        elif key.lower() == 'q':
            if self.selected_points:
                print("撤销上一次选择的点。")
                self.selected_points.pop()
                actor = self.point_actors.pop()
                self.renderer.RemoveActor(actor)
                self.render_window.Render()
            else:
                print("没有可撤销的点。")
        elif key.lower() == 'a':
            if self.centerlines:
                print("正在计算横截面...")
                self.calculate_cross_sections()
            else:
                print("请先生成中心线。")
        elif key.lower() == 's':
            if self.centerlines:
                print("正在计算表面面积...")
                self.calculate_surface_areas()
            else:
                print("请先生成中心线。")

    def display_point(self, position):
        # Create a sphere to represent the point
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(position)
        sphere.SetRadius(0.5)
        sphere.Update()

        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # Red color

        # Add the actor to the scene and store it
        self.renderer.AddActor(actor)
        self.point_actors.append(actor)
        self.render_window.Render()

    def generate_centerline(self):
        try:
            # Extract the source and target points as flat lists
            source_point = [float(coord) for coord in self.selected_points[0]]
            target_point = [float(coord) for coord in self.selected_points[1]]

            # 打印源点和目标点以检查
            print("源点：", source_point)
            print("目标点：", target_point)

            # Run vmtkcenterlines
            self.centerlines = vmtkscripts.vmtkCenterlines()
            self.centerlines.Surface = self.surface
            self.centerlines.SeedSelectorName = 'pointlist'
            self.centerlines.SourcePoints = source_point
            self.centerlines.TargetPoints = target_point
            self.centerlines.AppendEndPoints = 1
            self.centerlines.Resampling = 1
            self.centerlines.ResamplingStepLength = 0.5  # 可以根据需要调整步长
            self.centerlines.Execute()

            # 检查中心线是否生成成功
            if not self.centerlines.Centerlines or self.centerlines.Centerlines.GetNumberOfPoints() == 0:
                raise ValueError("中心线生成失败，请检查模型和选取的点。")

            # 对中心线进行平滑处理
            self.smooth_centerline()

            # 检查中心线的方向，必要时反转
            self.check_and_reverse_centerline()

            # Save the centerlines to a file
            writer = vtk.vtkPolyDataWriter()
            writer.SetFileName(self.output_file)
            writer.SetInputData(self.centerlines.Centerlines)
            writer.Write()
            print("中心线已保存到 {}".format(self.output_file))

            # Display the centerlines
            self.display_centerlines(self.centerlines.Centerlines)

            # 将表面模型设置为半透明
            self.surface_actor.GetProperty().SetOpacity(0.3)
            self.render_window.Render()
        except Exception as e:
            print("生成中心线时出错：{}".format(e))
            traceback.print_exc()
            # 不退出程序，允许用户继续操作

    def smooth_centerline(self):
        # 使用 vtkSplineFilter 对中心线进行平滑
        spline_filter = vtk.vtkSplineFilter()
        spline_filter.SetInputData(self.centerlines.Centerlines)
        spline_filter.SetSubdivideToLength()
        spline_filter.SetLength(0.5)  # 可以根据需要调整
        spline_filter.Update()
        self.centerlines.Centerlines = spline_filter.GetOutput()

    def check_and_reverse_centerline(self):
        centerline_points = self.centerlines.Centerlines.GetPoints()
        num_points = centerline_points.GetNumberOfPoints()
        first_centerline_point = centerline_points.GetPoint(0)
        last_centerline_point = centerline_points.GetPoint(num_points - 1)

        dist_to_first = np.linalg.norm(np.array(first_centerline_point) - np.array(self.selected_points[0]))
        dist_to_last = np.linalg.norm(np.array(last_centerline_point) - np.array(self.selected_points[0]))

        if dist_to_first > dist_to_last:
            # 反转中心线
            reversed_centerline = vtk.vtkPolyData()
            reversed_points = vtk.vtkPoints()
            reversed_lines = vtk.vtkCellArray()

            for i in range(num_points):
                reversed_points.InsertNextPoint(centerline_points.GetPoint(num_points - i - 1))

            reversed_centerline.SetPoints(reversed_points)

            line = vtk.vtkPolyLine()
            line.GetPointIds().SetNumberOfIds(num_points)
            for i in range(num_points):
                line.GetPointIds().SetId(i, i)

            reversed_lines.InsertNextCell(line)
            reversed_centerline.SetLines(reversed_lines)

            self.centerlines.Centerlines = reversed_centerline
            print("中心线已反转，以第一个选择的点为起点。")

    def display_centerlines(self, centerlines):
        # Remove existing centerline actor if present
        if self.centerline_actor:
            self.renderer.RemoveActor(self.centerline_actor)

        # Create a mapper and actor for the centerlines
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(centerlines)
        self.centerline_actor = vtk.vtkActor()
        self.centerline_actor.SetMapper(mapper)
        self.centerline_actor.GetProperty().SetColor(1.0, 1.0, 0.0)  # Yellow color
        self.centerline_actor.GetProperty().SetLineWidth(4)
        self.centerline_actor.GetProperty().SetOpacity(1.0)

        # Add the centerline actor to the scene
        self.renderer.AddActor(self.centerline_actor)
        self.render_window.Render()

    def calculate_cross_sections(self):
        if not self.centerlines:
            print("请先生成中心线。")
            return

        # Remove existing cross-section actors
        for actor in self.cross_section_actors:
            self.renderer.RemoveActor(actor)
        self.cross_section_actors = []
        self.cross_section_areas = []
        self.distance_along_curve = []
        self.cross_section_points = []

        # Get the centerline points
        points = self.centerlines.Centerlines.GetPoints()
        num_points = points.GetNumberOfPoints()

        if num_points < 3:
            print("中心线点数不足，无法计算横截面。")
            return

        # Calculate cumulative distance along the centerline
        cumulative_distances = [0.0]
        for i in range(1, num_points):
            p0 = np.array(points.GetPoint(i - 1))
            p1 = np.array(points.GetPoint(i))
            dist = np.linalg.norm(p1 - p0)
            cumulative_distances.append(cumulative_distances[-1] + dist)

        self.distance_along_curve = cumulative_distances

        # Preprocess the surface to remove non-manifold edges
        self.preprocess_surface()

        for i in range(num_points):
            point = np.array(points.GetPoint(i))
            self.cross_section_points.append(point)

            # Compute tangent vector using neighboring points
            if i == 0:
                p_prev = np.array(points.GetPoint(i))
                p_next = np.array(points.GetPoint(i + 1))
            elif i == num_points - 1:
                p_prev = np.array(points.GetPoint(i - 1))
                p_next = np.array(points.GetPoint(i))
            else:
                p_prev = np.array(points.GetPoint(i - 1))
                p_next = np.array(points.GetPoint(i + 1))

            tangent = p_next - p_prev

            # Normalize the tangent vector
            if np.linalg.norm(tangent) != 0:
                tangent = tangent / np.linalg.norm(tangent)
            else:
                tangent = np.array([1.0, 0.0, 0.0])

            # Create a plane perpendicular to the tangent at this point
            plane = vtk.vtkPlane()
            plane.SetOrigin(point)
            plane.SetNormal(tangent)

            # Perform boolean intersection to get the cross-section
            cross_section = self.get_cross_section(self.surface, plane, point, radius=5.0)

            # 检查横截面是否有效
            if cross_section.GetNumberOfPoints() == 0 or cross_section.GetNumberOfPolys() == 0:
                print("在索引 {} 处未能生成有效的横截面。".format(i))
                area = 0.0
                self.cross_section_areas.append(area)
                continue

            # 打印调试信息
            print("索引 {}：横截面有 {} 个点和 {} 个多边形。".format(
                i, cross_section.GetNumberOfPoints(), cross_section.GetNumberOfPolys()))

            # Calculate area
            area = self.calculate_area(cross_section)

            # Display the cross-section
            self.display_cross_section(cross_section)

            self.cross_section_areas.append(area)

        print("横截面计算完成。")

    def get_cross_section(self, surface, plane, point, radius=5.0):
        # 创建一个球形裁剪器，以限制横截面的范围
        sphere = vtk.vtkSphere()
        sphere.SetCenter(point)
        sphere.SetRadius(radius)

        # 使用 vtkClipPolyData 来裁剪表面，保留球体内部的部分
        clipper_sphere = vtk.vtkClipPolyData()
        clipper_sphere.SetInputData(surface)
        clipper_sphere.SetClipFunction(sphere)
        clipper_sphere.InsideOutOn()
        clipper_sphere.Update()
        clipped_surface = clipper_sphere.GetOutput()

        # 使用裁剪后的表面进行截面计算
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputData(clipped_surface)
        cutter.Update()

        cross_section = cutter.GetOutput()

        # 检查是否有足够的点来形成多边形
        if cross_section.GetNumberOfPoints() < 3:
            return vtk.vtkPolyData()

        # 使用 vtkStripper 来连接线段，形成闭合的多边形
        stripper = vtk.vtkStripper()
        stripper.SetInputData(cross_section)
        stripper.JoinContiguousSegmentsOn()
        stripper.Update()

        # 检查是否有足够的点
        if stripper.GetOutput().GetNumberOfPoints() < 3:
            return vtk.vtkPolyData()

        # 使用 vtkContourTriangulator 来生成多边形
        contour = vtk.vtkContourTriangulator()
        contour.SetInputData(stripper.GetOutput())
        contour.Update()

        triangulated = contour.GetOutput()

        return triangulated

    def calculate_area(self, polydata):
        if polydata.GetNumberOfPolys() == 0:
            return 0.0
        mass = vtk.vtkMassProperties()
        mass.SetInputData(polydata)
        area = mass.GetSurfaceArea()
        return area

    def display_cross_section(self, cross_section):
        # Create a mapper and actor for the cross-section
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(cross_section)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1.0, 0.0, 1.0)  # Magenta color
        actor.GetProperty().SetLineWidth(2)
        actor.GetProperty().SetOpacity(1.0)
        actor.GetProperty().SetRepresentationToSurface()

        # Add the actor to the renderer
        self.renderer.AddActor(actor)
        self.cross_section_actors.append(actor)
        self.render_window.Render()

    def calculate_surface_areas(self):
        if not self.cross_section_areas:
            print("请先计算横截面。")
            return

        num_points = len(self.cross_section_points)
        self.surface_areas = []
        for i in range(num_points):
            if i == 0:
                area = 0.0
            else:
                # 计算两个相邻横截面之间的表面积段
                distance = self.distance_along_curve[i] - self.distance_along_curve[i - 1]
                avg_perimeter = (self.cross_section_areas[i] + self.cross_section_areas[i - 1]) / 2.0
                area = avg_perimeter * distance
            self.surface_areas.append(area)
        print("表面面积计算完成。")

    def save_to_csv(self, filename):
        # Save the collected data to a CSV file
        headers = ['distance along curve', 'x', 'y', 'z', 'cross-sectional area at vertex', 'surface area between points']
        with open(filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(headers)
            num_points = len(self.cross_section_points)
            for i in range(num_points):
                point = self.cross_section_points[i]
                distance = self.distance_along_curve[i]
                area = self.cross_section_areas[i] if i < len(self.cross_section_areas) else 0.0
                surface_area = self.surface_areas[i] if i < len(self.surface_areas) else 0.0
                csvwriter.writerow([distance, point[0], point[1], point[2], area, surface_area])

        print("数据已保存到 {}".format(filename))

    def run(self):
        self.load_surface()
        self.setup_render()
        self.render_window.Render()
        self.interactor.Start()

        # After interactor ends, check if data is available
        if self.cross_section_points and self.distance_along_curve:
            # Save data to CSV
            csv_filename = os.path.splitext(self.output_file)[0] + '_data.csv'
            self.save_to_csv(csv_filename)
        else:
            print("没有数据需要保存。")

if __name__ == '__main__':
    # 请将以下文件路径替换为您的实际文件路径
    surface_file = "C:/Users/qd261/Desktop/PD_Study/Secretin_MRCP_Simple/3-Matic_Model/6 NS.stl"
    output_file = "C:/Users/qd261/Desktop/6 NS_centerline.dat"

    centerline_extractor = CenterlineExtraction(surface_file, output_file)
    centerline_extractor.run()
