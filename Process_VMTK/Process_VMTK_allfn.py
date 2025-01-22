import vtk
from vmtk import vmtkscripts
import csv
import os
import sys
import traceback
import numpy as np

class CenterlineExtraction:
    def __init__(self, surface_file, output_file, **kwargs):
        self.surface_file = surface_file
        self.output_file = output_file
        self.surface = None
        self.renderer = None
        self.render_window = None
        self.interactor = None
        self.selected_points = []
        self.point_actors = []
        self.picker = vtk.vtkCellPicker()
        # 使用参数中的picker_tolerance，默认为0.005
        self.picker_tolerance = kwargs.get('picker_tolerance', 0.005)
        self.picker.SetTolerance(self.picker_tolerance)
        self.centerline_actor = None
        self.surface_actor = None
        self.centerlines = None
        self.cross_section_actors = []
        self.cross_section_areas = []
        self.distance_along_curve = []
        self.cross_section_points = []
        self.surface_areas = []
        # 其他可调参数
        self.resampling_step_length = kwargs.get('resampling_step_length', 0.1)
        self.spline_filter_length = kwargs.get('spline_filter_length', 0.1)
        self.cross_section_radius = kwargs.get('cross_section_radius', 15.0)
        self.seed_selector_name = kwargs.get('seed_selector_name', 'pointlist')
        self.append_end_points = kwargs.get('append_end_points', 1)
        self.resampling = kwargs.get('resampling', 1)
        self.sphere_radius = kwargs.get('sphere_radius', 5.0)
        self.sphere_inside_out = kwargs.get('sphere_inside_out', True)
        self.cross_section_display_color = kwargs.get('cross_section_display_color', (1.0, 0.0, 1.0))
        self.cross_section_line_width = kwargs.get('cross_section_line_width', 2)

    def load_surface(self):
        try:
            # 读取STL文件
            reader = vtk.vtkSTLReader()
            reader.SetFileName(self.surface_file)
            reader.Update()
            self.surface = reader.GetOutput()
            if self.surface.GetNumberOfPoints() == 0:
                raise ValueError("The surface model could not be loaded correctly. Please check the file path and file format.")
        except Exception as e:
            print("Error loading surface model: {}".format(e))
            traceback.print_exc()

    def preprocess_surface(self):
        # 使用vtkCleanPolyData清理表面，移除非流形几何
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(self.surface)
        cleaner.Update()
        self.surface = cleaner.GetOutput()

    def setup_render(self):
        # 创建表面的mapper和actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.surface)
        self.surface_actor = vtk.vtkActor()
        self.surface_actor.SetMapper(mapper)

        # 创建渲染器、渲染窗口和交互器
        self.renderer = vtk.vtkRenderer()
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.render_window)

        # 将表面actor添加到场景中
        self.renderer.AddActor(self.surface_actor)
        self.renderer.SetBackground(0.1, 0.2, 0.4)

        # 设置picker
        self.interactor.SetPicker(self.picker)

        # 使用旋转的交互风格
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        # 设置事件处理
        self.interactor.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.interactor.AddObserver("KeyPressEvent", self.key_press_event)

    def left_button_press_event(self, obj, event):
        click_pos = self.interactor.GetEventPosition()
        self.picker.Pick(click_pos[0], click_pos[1], 0, self.renderer)
        world_position = self.picker.GetPickPosition()
        cell_id = self.picker.GetCellId()
        if cell_id >= 0:
            # 将选择的点添加到列表中
            self.selected_points.append(world_position)
            print("Selected point coordinates: {}".format(world_position))
            # 显示选择的点
            self.display_point(world_position)
        else:
            print("No valid surface point was selected.")

    def key_press_event(self, obj, event):
        key = self.interactor.GetKeySym()
        if key == 'space':
            if len(self.selected_points) == 2:
                print("Two points selected, generating centerline...")
                self.generate_centerline()
            else:
                print("Selected {} points. Please select two points to continue.".format(len(self.selected_points)))
        elif key.lower() == 'q':
            if self.selected_points:
                print("Undoing the last selected point.")
                self.selected_points.pop()
                actor = self.point_actors.pop()
                self.renderer.RemoveActor(actor)
                self.render_window.Render()
            else:
                print("No points to undo.")
        elif key.lower() == 'a':
            if self.centerlines:
                print("Calculating cross sections...")
                self.calculate_cross_sections()
            else:
                print("Please generate the centerline first.")
        elif key.lower() == 's':
            if self.centerlines:
                print("Calculating surface areas...")
                self.calculate_surface_areas()
            else:
                print("Please generate the centerline first.")

    def display_point(self, position):
        # 创建一个球体来表示点
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(position)
        sphere.SetRadius(0.5)
        sphere.Update()

        # 创建mapper和actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # 红色

        # 将actor添加到场景中并存储
        self.renderer.AddActor(actor)
        self.point_actors.append(actor)
        self.render_window.Render()

    def generate_centerline(self):
        try:
            # 提取源点和目标点
            source_point = [float(coord) for coord in self.selected_points[0]]
            target_point = [float(coord) for coord in self.selected_points[1]]

            # 打印源点和目标点以验证
            print("Source point:", source_point)
            print("Target point:", target_point)

            # 运行vmtkcenterlines
            self.centerlines = vmtkscripts.vmtkCenterlines()
            self.centerlines.Surface = self.surface
            self.centerlines.SeedSelectorName = self.seed_selector_name
            self.centerlines.SourcePoints = source_point
            self.centerlines.TargetPoints = target_point
            self.centerlines.AppendEndPoints = self.append_end_points
            self.centerlines.Resampling = self.resampling
            self.centerlines.ResamplingStepLength = self.resampling_step_length  # 使用参数
            self.centerlines.Execute()

            # 检查中心线是否成功生成
            if not self.centerlines.Centerlines or self.centerlines.Centerlines.GetNumberOfPoints() == 0:
                raise ValueError("Centerline generation failed. Please check the model and selected points.")

            # 平滑中心线
            self.smooth_centerline()

            # 检查并在必要时反转中心线方向
            self.check_and_reverse_centerline()

            # 将中心线保存到文件
            writer = vtk.vtkPolyDataWriter()
            writer.SetFileName(self.output_file)
            writer.SetInputData(self.centerlines.Centerlines)
            writer.Write()
            print("Centerline has been saved to {}".format(self.output_file))

            # 显示中心线
            self.display_centerlines(self.centerlines.Centerlines)

            # 将表面模型设置为半透明
            self.surface_actor.GetProperty().SetOpacity(0.3)
            self.render_window.Render()
        except Exception as e:
            print("Error generating centerline: {}".format(e))
            traceback.print_exc()

    def smooth_centerline(self):
        # 使用vtkSplineFilter平滑中心线
        spline_filter = vtk.vtkSplineFilter()
        spline_filter.SetInputData(self.centerlines.Centerlines)
        spline_filter.SetSubdivideToLength()
        spline_filter.SetLength(self.spline_filter_length)  # 使用参数
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
            print("Centerline has been reversed to start from the first selected point.")

    def display_centerlines(self, centerlines):
        # 如果存在现有的中心线actor，移除它
        if self.centerline_actor:
            self.renderer.RemoveActor(self.centerline_actor)

        # 创建中心线的mapper和actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(centerlines)
        self.centerline_actor = vtk.vtkActor()
        self.centerline_actor.SetMapper(mapper)
        self.centerline_actor.GetProperty().SetColor(1.0, 1.0, 0.0)  # 黄色
        self.centerline_actor.GetProperty().SetLineWidth(4)
        self.centerline_actor.GetProperty().SetOpacity(1.0)

        # 将中心线actor添加到场景中
        self.renderer.AddActor(self.centerline_actor)
        self.render_window.Render()

    def calculate_cross_sections(self):
        if not self.centerlines:
            print("Please generate the centerline first.")
            return

        # 移除现有的截面actor
        for actor in self.cross_section_actors:
            self.renderer.RemoveActor(actor)
        self.cross_section_actors = []
        self.cross_section_areas = []
        self.distance_along_curve = []
        self.cross_section_points = []

        # 获取中心线点
        points = self.centerlines.Centerlines.GetPoints()
        num_points = points.GetNumberOfPoints()

        if num_points < 3:
            print("Insufficient centerline points to calculate cross sections.")
            return

        # 计算沿中心线的累积距离
        cumulative_distances = [0.0]
        for i in range(1, num_points):
            p0 = np.array(points.GetPoint(i - 1))
            p1 = np.array(points.GetPoint(i))
            dist = np.linalg.norm(p1 - p0)
            cumulative_distances.append(cumulative_distances[-1] + dist)

        self.distance_along_curve = cumulative_distances

        # 预处理表面以移除非流形边
        self.preprocess_surface()

        for i in range(num_points):
            point = np.array(points.GetPoint(i))
            self.cross_section_points.append(point)

            # 使用邻近点计算切向量
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

            # 归一化切向量
            if np.linalg.norm(tangent) != 0:
                tangent = tangent / np.linalg.norm(tangent)
            else:
                tangent = np.array([1.0, 0.0, 0.0])

            # 创建在该点垂直于切向量的平面
            plane = vtk.vtkPlane()
            plane.SetOrigin(point)
            plane.SetNormal(tangent)

            # 执行布尔交集以获取截面
            cross_section = self.get_cross_section(self.surface, plane, point)

            # 检查截面是否有效
            if cross_section.GetNumberOfPoints() == 0 or cross_section.GetNumberOfPolys() == 0:
                print("Failed to generate a valid cross section at index {}.".format(i))
                area = 0.0
                self.cross_section_areas.append(area)
                continue

            # 打印调试信息
            print("Index {}: Cross section has {} points and {} polygons.".format(
                i, cross_section.GetNumberOfPoints(), cross_section.GetNumberOfPolys()))

            # 计算面积
            area = self.calculate_area(cross_section)

            # 显示截面
            self.display_cross_section(cross_section)

            self.cross_section_areas.append(area)

        print("Cross section calculation completed.")

    def get_cross_section(self, surface, plane, point):
        # 创建一个球形裁剪器以限制截面范围
        sphere = vtk.vtkSphere()
        sphere.SetCenter(point)
        sphere.SetRadius(self.cross_section_radius)

        # 使用vtkClipPolyData裁剪表面，保留球体内部
        clipper_sphere = vtk.vtkClipPolyData()
        clipper_sphere.SetInputData(surface)
        clipper_sphere.SetClipFunction(sphere)
        if self.sphere_inside_out:
            clipper_sphere.InsideOutOn()
        else:
            clipper_sphere.InsideOutOff()
        clipper_sphere.Update()
        clipped_surface = clipper_sphere.GetOutput()

        # 使用裁剪后的表面进行截面计算
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputData(clipped_surface)
        cutter.Update()

        cross_section = cutter.GetOutput()

        # 检查是否有足够的点形成多边形
        if cross_section.GetNumberOfPoints() < 3:
            return vtk.vtkPolyData()

        # 使用vtkStripper将线段连接成闭合多边形
        stripper = vtk.vtkStripper()
        stripper.SetInputData(cross_section)
        stripper.JoinContiguousSegmentsOn()
        stripper.Update()

        # 检查是否有足够的点
        if stripper.GetOutput().GetNumberOfPoints() < 3:
            return vtk.vtkPolyData()

        # 使用vtkContourTriangulator生成多边形
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
        # 创建截面的mapper和actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(cross_section)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*self.cross_section_display_color)  # 使用参数
        actor.GetProperty().SetLineWidth(self.cross_section_line_width)
        actor.GetProperty().SetOpacity(1.0)
        actor.GetProperty().SetRepresentationToSurface()

        # 将actor添加到渲染器中
        self.renderer.AddActor(actor)
        self.cross_section_actors.append(actor)
        self.render_window.Render()

    def calculate_surface_areas(self):
        if not self.cross_section_areas:
            print("Please calculate cross sections first.")
            return

        num_points = len(self.cross_section_points)
        self.surface_areas = []
        for i in range(num_points):
            if i == 0:
                area = 0.0
            else:
                # 计算两个相邻截面之间的表面积段
                distance = self.distance_along_curve[i] - self.distance_along_curve[i - 1]
                avg_perimeter = (self.cross_section_areas[i] + self.cross_section_areas[i - 1]) / 2.0
                area = avg_perimeter * distance
            self.surface_areas.append(area)
        print("Surface area calculation completed.")

    def save_to_csv(self, filename):
        # 将收集的数据保存到CSV文件
        headers = ['distance along curve', 'x', 'y', 'z', 'cross-sectional area at vertex', 'sum of surface areas projected to this vertex']
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

        print("Data has been saved to {}".format(filename))

    def run(self):
        self.load_surface()
        self.setup_render()
        self.render_window.Render()
        self.interactor.Start()

        # 在交互结束后，检查数据是否可用
        if self.cross_section_points and self.distance_along_curve:
            # 将数据保存到CSV
            csv_filename = os.path.splitext(self.output_file)[0] + '.csv'
            self.save_to_csv(csv_filename)
        else:
            print("No data to save.")

if __name__ == '__main__':
    # 请将以下文件路径替换为您的实际文件路径

    N = 8
    dict = {'1': '1 DJ', 
            '2':'2 MD',  
            '4':'4 SM', 
            '6':'6 NS',
            '8': '8 PM', 
            '13':'13 LB', 
            '14':'14 VK', 
            '15':'15 SK', 
            '16':'16 GGG', 
            '17': '17 VJ', 
            '18':'18 MDS', 
            '19':'19 CJ', 
            '20':'20 TA', 
            '21':'21 MSA', 
            '24':'24 BKS',
            '26':'26 NR', 
            '27':'27 NA'}

    # ## Secretin_MRCP_GROUP
    # N = input("file #:")
    # input_name = dict[str(N)]
    # output_name = input_name
    # surface_file = "C:/Users/qd261/Desktop/PD_Study/Secretin_MRCP_Simple/3-Matic_Model/{}.stl".format(input_name)
    # output_file = r'C:\Users\qd261\Desktop\PD_Study\Secretin_MRCP_Simple\VMTK\{}.scv'.format(output_name)



    ## Secretin_MRCP_GROUP
    # N = input("file #:")
    # input_name = dict[str(N)]
    # output_name = input_name



    # surface_file = r"C:\Users\qd261\Desktop\Secretin_MRCP_Simple_new\3-Matic\08_Exp84.stl"
    # output_file = r'C:\Users\qd261\Desktop\08_Exp84.csv'


    Name = '27'
    surface_file = r"C:\Users\qd261\Desktop\Secretin_MRCP_Simple_new\3-Matic\{}.stl".format(Name)
    output_file = r'C:\Users\qd261\Desktop\Secretin_MRCP_Simple_new\centerline\{}.csv'.format(Name)

    # ## Hopkins CP1
    # N = input("file #:")
    # input_name = str(N)
    # output_name = input_name
    # surface_file =r"C:\Users\qd261\Desktop\Hopkins CP1-REDO\3-matic\{}.stl".format(input_name)
    # output_file = r'C:\Users\qd261\Desktop\Hopkins CP1-REDO\VMTK\{}.scv'.format(output_name)

    # ## Hopkins CP2
    # N = input("file #:")
    # input_name = str(N)
    # output_name = input_name
    # surface_file =r"C:\Users\qd261\Desktop\Hopkins CP2-REDO\3-Matic\{}.stl".format(input_name)
    # output_file = r'C:\Users\qd261\Desktop\Hopkins CP2-REDO\VMTK\{}.scv'.format(output_name)


    
    
    # 设置可调参数
    params = {
        'picker_tolerance': 0.005,
        'resampling_step_length': 0.05,
        'spline_filter_length': 0.5,
        'cross_section_radius': 20.0,
        'seed_selector_name': 'pointlist',
        'append_end_points': 1,
        'resampling': 1,
        'sphere_radius': 5.0,
        'sphere_inside_out': True,
        'cross_section_display_color': (1.0, 0.0, 1.0),
        'cross_section_line_width': 2,
    }

    centerline_extractor = CenterlineExtraction(surface_file, output_file, **params)
    centerline_extractor.run()
