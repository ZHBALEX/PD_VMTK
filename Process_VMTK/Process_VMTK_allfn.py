import vtk
from vmtk import vmtkscripts
import csv
import os
import sys
import traceback

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
        self.previous_valid_tangent = None  # 保存前一个有效的切向量

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
            sys.exit(1)

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

            # 不进行中心线平滑，直接使用生成的中心线

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
            sys.exit(1)

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

        if num_points < 2:
            print("中心线点数不足，无法计算横截面。")
            return

        # Loop through each point on the centerline
        for i in range(num_points):
            point = points.GetPoint(i)
            self.cross_section_points.append(point)

            # Calculate the normal (tangent) at this point
            tangent = self.calculate_tangent(i, points)
            if tangent is None:
                # 使用前一个有效的切向量
                if self.previous_valid_tangent is not None:
                    tangent = self.previous_valid_tangent
                    print("使用前一个有效的切向量，索引 {}".format(i))
                else:
                    # 无法计算切向量，跳过该点
                    print("无法计算切向量，跳过索引 {}".format(i))
                    self.cross_section_areas.append(0.0)
                    if i == 0:
                        distance = 0.0
                    else:
                        distance = self.distance_along_curve[-1]
                    self.distance_along_curve.append(distance)
                    continue
            else:
                # 保存当前有效的切向量
                self.previous_valid_tangent = tangent

            # Create a plane perpendicular to the tangent at this point
            plane = vtk.vtkPlane()
            plane.SetOrigin(point)
            plane.SetNormal(tangent)

            # Clip the surface with the plane to get the cross-section
            cutter = vtk.vtkCutter()
            cutter.SetCutFunction(plane)
            cutter.SetInputData(self.surface)
            cutter.GenerateValues(1, 0, 0)
            cutter.Update()

            cross_section = cutter.GetOutput()

            # Ensure cross_section is valid
            if cross_section.GetNumberOfPoints() == 0:
                print("在索引 {} 处未能生成有效的横截面。".format(i))
                area = 0.0
                self.cross_section_areas.append(area)
                # Calculate distance along the curve
                if i == 0:
                    distance = 0.0
                else:
                    prev_point = points.GetPoint(i - 1)
                    distance = self.distance_along_curve[-1] + self.euclidean_distance(point, prev_point)
                self.distance_along_curve.append(distance)
                continue

            # Process the cross-section to create a polygon and triangulate it
            # Step 1: Use vtkStripper to create a closed polyline
            stripper = vtk.vtkStripper()
            stripper.SetInputData(cross_section)
            stripper.Update()

            # Step 2: Create a polygon from the polyline
            contour = vtk.vtkContourTriangulator()
            contour.SetInputData(stripper.GetOutput())
            contour.Update()

            triangulated = contour.GetOutput()

            # Check if triangulation was successful
            if triangulated.GetNumberOfCells() == 0:
                print("在索引 {} 处无法三角化横截面。".format(i))
                area = 0.0
                self.cross_section_areas.append(area)
                # Calculate distance along the curve
                if i == 0:
                    distance = 0.0
                else:
                    prev_point = points.GetPoint(i - 1)
                    distance = self.distance_along_curve[-1] + self.euclidean_distance(point, prev_point)
                self.distance_along_curve.append(distance)
                continue

            area = self.calculate_area(triangulated)

            # Display the cross-section
            self.display_cross_section(triangulated)

            self.cross_section_areas.append(area)

            # Calculate distance along the curve
            if i == 0:
                distance = 0.0
            else:
                prev_point = points.GetPoint(i - 1)
                distance = self.distance_along_curve[-1] + self.euclidean_distance(point, prev_point)
            self.distance_along_curve.append(distance)

        print("横截面计算完成。")

    def calculate_tangent(self, index, points):
        num_points = points.GetNumberOfPoints()
        if num_points < 2:
            return None

        # 使用加权平均方法计算切向量，考虑前后各两个点
        p0 = points.GetPoint(index)

        if index == 0 or index == 1:
            p_m1 = points.GetPoint(1)
            p_m2 = points.GetPoint(0)
            norm_m1 = [p_m1[i] - p_m2[i] for i in range(3)]
            norm_m2 = norm_m1
        else:
            p_m1 = points.GetPoint(index - 1)
            p_m2 = points.GetPoint(index - 2)
            norm_m1 = [p0[i] - p_m1[i] for i in range(3)]
            norm_m2 = [p_m1[i] - p_m2[i] for i in range(3)]

        if index == num_points - 1 or index == num_points - 2:
            p_p1 = points.GetPoint(num_points - 1)
            p_p2 = points.GetPoint(num_points - 2)
            norm_p1 = [p_p1[i] - p_p2[i] for i in range(3)]
            norm_p2 = norm_p1
        else:
            p_p1 = points.GetPoint(index + 1)
            p_p2 = points.GetPoint(index + 2)
            norm_p1 = [p_p1[i] - p0[i] for i in range(3)]
            norm_p2 = [p_p2[i] - p_p1[i] for i in range(3)]

        # 计算加权平均
        norm_m1 = self.normalize_vector(norm_m1)
        norm_m2 = self.normalize_vector(norm_m2)
        norm_p1 = self.normalize_vector(norm_p1)
        norm_p2 = self.normalize_vector(norm_p2)

        tangent = [norm_m1[i] + norm_p1[i] + 0.5 * (norm_m2[i] + norm_p2[i]) for i in range(3)]

        norm = sum([tangent[i] ** 2 for i in range(3)]) ** 0.5
        if norm == 0:
            print("切向量为零，跳过索引 {}".format(index))
            return None
        tangent = [tangent[i] / norm for i in range(3)]
        return tangent

    def normalize_vector(self, vec):
        norm = sum([vec[i] ** 2 for i in range(3)]) ** 0.5
        if norm == 0:
            return [0.0, 0.0, 0.0]
        return [vec[i] / norm for i in range(3)]

    def calculate_area(self, polydata):
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
                distance = self.euclidean_distance(self.cross_section_points[i], self.cross_section_points[i - 1])
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
                if i < len(self.cross_section_points):
                    point = self.cross_section_points[i]
                else:
                    point = (0.0, 0.0, 0.0)
                if i < len(self.distance_along_curve):
                    distance = self.distance_along_curve[i]
                else:
                    distance = 0.0
                if i < len(self.cross_section_areas):
                    area = self.cross_section_areas[i]
                else:
                    area = 0.0
                if i < len(self.surface_areas):
                    surface_area = self.surface_areas[i]
                else:
                    surface_area = 0.0
                csvwriter.writerow([distance, point[0], point[1], point[2], area, surface_area])

        print("数据已保存到 {}".format(filename))

    def euclidean_distance(self, p1, p2):
        return sum([(p1[i] - p2[i]) ** 2 for i in range(3)]) ** 0.5

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
