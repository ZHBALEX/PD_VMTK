import vtk
from vmtk import vmtkscripts

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
        self.picker.SetTolerance(0.0005)
        self.centerline_actor = None  # 用于存储中心线的actor
        self.surface_actor = None     # 用于存储表面模型的actor

    def load_surface(self):
        # 读取STL文件
        reader = vtk.vtkSTLReader()
        reader.SetFileName(self.surface_file)
        reader.Update()
        self.surface = reader.GetOutput()

    def setup_render(self):
        # 为表面创建mapper和actor
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

        # 使用TrackballCamera风格进行旋转
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
            # 添加选取的点到列表
            self.selected_points.append(world_position)
            print("选择的点坐标: {}".format(world_position))
            # 显示选取的点
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

    def display_point(self, position):
        # 创建球体来表示选取的点
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
        # 提取源点和目标点为浮点数列表
        source_point = [float(coord) for coord in self.selected_points[0]]
        target_point = [float(coord) for coord in self.selected_points[1]]

        # 运行vmtkcenterlines
        centerlines = vmtkscripts.vmtkCenterlines()
        centerlines.Surface = self.surface
        centerlines.SeedSelectorName = 'pointlist'
        centerlines.SourcePoints = source_point
        centerlines.TargetPoints = target_point
        centerlines.AppendEndPoints = 1
        centerlines.Resampling = 1
        centerlines.ResamplingStepLength = 0.1
        centerlines.Execute()

        # 将中心线保存到文件
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(self.output_file)
        writer.SetInputData(centerlines.Centerlines)
        writer.Write()
        print("中心线已保存到 {}".format(self.output_file))

        # 显示中心线
        self.display_centerlines(centerlines.Centerlines)

        # 将表面模型设置为半透明
        self.surface_actor.GetProperty().SetOpacity(0.3)
        self.render_window.Render()

    def display_centerlines(self, centerlines):
        # 如果存在已显示的中心线actor，先移除
        if self.centerline_actor:
            self.renderer.RemoveActor(self.centerline_actor)

        # 创建mapper和actor用于中心线
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

    def run(self):
        self.load_surface()
        self.setup_render()
        self.render_window.Render()
        self.interactor.Start()

if __name__ == '__main__':
    # 替换为你的文件路径
    surface_file = "C:/Users/qd261/Desktop/PD_Study/Secretin_MRCP_Simple/3-Matic_Model/6 NS.stl"
    output_file = "C:/Users/qd261/Desktop/6 NS_centerline.dat"

    centerline_extractor = CenterlineExtraction(surface_file, output_file)
    centerline_extractor.run()
