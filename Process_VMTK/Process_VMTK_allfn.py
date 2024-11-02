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
        self.picker.SetTolerance(0.005)  # Increase tolerance value
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
                raise ValueError("The surface model could not be loaded correctly. Please check the file path and file format.")
        except Exception as e:
            print("Error loading surface model: {}".format(e))
            traceback.print_exc()
            # Do not exit the program, allow the user to continue

    def preprocess_surface(self):
        # Use vtkCleanPolyData to clean the surface, removing non-manifold geometry
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
            print("Selected point coordinates: {}".format(world_position))
            # Display the picked point
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

            # Print source and target points for verification
            print("Source point:", source_point)
            print("Target point:", target_point)

            # Run vmtkcenterlines
            self.centerlines = vmtkscripts.vmtkCenterlines()
            self.centerlines.Surface = self.surface
            self.centerlines.SeedSelectorName = 'pointlist'
            self.centerlines.SourcePoints = source_point
            self.centerlines.TargetPoints = target_point
            self.centerlines.AppendEndPoints = 1
            self.centerlines.Resampling = 1
            self.centerlines.ResamplingStepLength = 0.1  # Adjust the step length as needed
            self.centerlines.Execute()

            # Check if centerlines were successfully generated
            if not self.centerlines.Centerlines or self.centerlines.Centerlines.GetNumberOfPoints() == 0:
                raise ValueError("Centerline generation failed. Please check the model and selected points.")

            # Smooth the centerline
            self.smooth_centerline()

            # Check and reverse the direction of the centerline if necessary
            self.check_and_reverse_centerline()

            # Save the centerlines to a file
            writer = vtk.vtkPolyDataWriter()
            writer.SetFileName(self.output_file)
            writer.SetInputData(self.centerlines.Centerlines)
            writer.Write()
            print("Centerline has been saved to {}".format(self.output_file))

            # Display the centerlines
            self.display_centerlines(self.centerlines.Centerlines)

            # Set the surface model to semi-transparent
            self.surface_actor.GetProperty().SetOpacity(0.3)
            self.render_window.Render()
        except Exception as e:
            print("Error generating centerline: {}".format(e))
            traceback.print_exc()
            # Do not exit the program, allow the user to continue

    def smooth_centerline(self):
        # Use vtkSplineFilter to smooth the centerline
        spline_filter = vtk.vtkSplineFilter()
        spline_filter.SetInputData(self.centerlines.Centerlines)
        spline_filter.SetSubdivideToLength()
        spline_filter.SetLength(0.1)  # Adjust as needed
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
            # Reverse the centerline
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
            print("Please generate the centerline first.")
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
            print("Insufficient centerline points to calculate cross sections.")
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
            cross_section = self.get_cross_section(self.surface, plane, point, radius=15.0)

            # Check if the cross-section is valid
            if cross_section.GetNumberOfPoints() == 0 or cross_section.GetNumberOfPolys() == 0:
                print("Failed to generate a valid cross section at index {}.".format(i))
                area = 0.0
                self.cross_section_areas.append(area)
                continue

            # Print debugging information
            print("Index {}: Cross section has {} points and {} polygons.".format(
                i, cross_section.GetNumberOfPoints(), cross_section.GetNumberOfPolys()))

            # Calculate area
            area = self.calculate_area(cross_section)

            # Display the cross-section
            self.display_cross_section(cross_section)

            self.cross_section_areas.append(area)

        print("Cross section calculation completed.")

    def get_cross_section(self, surface, plane, point, radius=5.0):
        # Create a spherical clipper to limit the cross-section range
        sphere = vtk.vtkSphere()
        sphere.SetCenter(point)
        sphere.SetRadius(radius)

        # Use vtkClipPolyData to clip the surface, keeping the inside of the sphere
        clipper_sphere = vtk.vtkClipPolyData()
        clipper_sphere.SetInputData(surface)
        clipper_sphere.SetClipFunction(sphere)
        clipper_sphere.InsideOutOn()
        clipper_sphere.Update()
        clipped_surface = clipper_sphere.GetOutput()

        # Use the clipped surface for cross-section calculation
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputData(clipped_surface)
        cutter.Update()

        cross_section = cutter.GetOutput()

        # Check if there are enough points to form a polygon
        if cross_section.GetNumberOfPoints() < 3:
            return vtk.vtkPolyData()

        # Use vtkStripper to connect line segments into a closed polygon
        stripper = vtk.vtkStripper()
        stripper.SetInputData(cross_section)
        stripper.JoinContiguousSegmentsOn()
        stripper.Update()

        # Check if there are enough points
        if stripper.GetOutput().GetNumberOfPoints() < 3:
            return vtk.vtkPolyData()

        # Use vtkContourTriangulator to generate polygons
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
            print("Please calculate cross sections first.")
            return

        num_points = len(self.cross_section_points)
        self.surface_areas = []
        for i in range(num_points):
            if i == 0:
                area = 0.0
            else:
                # Calculate the surface area segment between two adjacent cross sections
                distance = self.distance_along_curve[i] - self.distance_along_curve[i - 1]
                avg_perimeter = (self.cross_section_areas[i] + self.cross_section_areas[i - 1]) / 2.0
                area = avg_perimeter * distance
            self.surface_areas.append(area)
        print("Surface area calculation completed.")

    def save_to_csv(self, filename):
        # Save the collected data to a CSV file
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

        # After interactor ends, check if data is available
        if self.cross_section_points and self.distance_along_curve:
            # Save data to CSV
            csv_filename = os.path.splitext(self.output_file)[0] + '.csv'
            self.save_to_csv(csv_filename)
        else:
            print("No data to save.")

if __name__ == '__main__':
    # Please replace the following file paths with your actual file paths

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
    N = input("file #:")
    input_name = dict[str(N)]
    output_name = input_name
    surface_file = "C:/Users/qd261/Desktop/PD_Study/Secretin_MRCP_Simple/3-Matic_Model/{}.stl".format(input_name)
    output_file = r'C:\Users\qd261\Desktop\PD_Study\Secretin_MRCP_Simple\VMTK\{}.scv'.format(output_name)

    output_file = r'C:\Users\qd261\Desktop\{}.scv'.format(output_name)
    centerline_extractor = CenterlineExtraction(surface_file, output_file)
    centerline_extractor.run()
