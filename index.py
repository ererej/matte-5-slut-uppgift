from manim import *
import numpy as np
angle = PI # the angle taken out of the semicircle
radius = 1

class TestScene(ThreeDScene):
    def construct(self):
        
        # angle = PI # the angle taken out of the semicircle
        # radius = 1


        temp_circle = Circle(radius=1, color=BLUE, fill_opacity=0.5)
        self.play(Create(temp_circle))
        


        upper_semicircleArc = Arc(radius=1, start_angle=0, angle=angle)
        lower_semicircleArc = Arc(radius=1, start_angle=angle, angle=2*PI - angle)

        upper_semicircle = VGroup(upper_semicircleArc)
        lower_semicircle = VGroup(lower_semicircleArc)

        # Create a group with the 2 semicircles
        semi_circles = VGroup(upper_semicircleArc, lower_semicircleArc)
        # Set the color of the semicircles
        semi_circles.set_color(BLUE)


        # Add the group to the scene
        self.add(upper_semicircleArc)
        self.add(lower_semicircleArc)
        self.wait(1)
        # split the lines into 2 to allow for other angles then 180 degrees
        upper_seperator_1 = Line(start=upper_semicircleArc.get_start(), end=upper_semicircleArc.get_arc_center())
        upper_seperator_2 = Line(start=upper_semicircleArc.get_start(), end=upper_semicircleArc.get_arc_center())

        lower_seperator_1 = Line(start=lower_semicircleArc.get_start(), end=lower_semicircleArc.get_arc_center())
        lower_seperator_2 = Line(start=lower_semicircleArc.get_arc_center(), end=lower_semicircleArc.get_end())  
        
        upper_seperator = VGroup(upper_seperator_1, upper_seperator_2)
        lower_seperator = VGroup(lower_seperator_1, lower_seperator_2)

        upper_semicircle.add(upper_seperator_1, upper_seperator_2, upper_semicircleArc)
        lower_semicircle.add(lower_seperator_1, lower_seperator_2, lower_semicircleArc)

        # Fill the upper semicircle
        upper_fill = AnnularSector(
            inner_radius=0,
            outer_radius=1,
            angle=angle,
            start_angle=0,
            color=BLUE,
            fill_opacity=0.5
        )
        self.add(upper_fill)
        upper_semicircle.add(upper_fill)

        # Fill the lower semicircle
        lower_fill = AnnularSector(
            inner_radius=0,
            outer_radius=1,
            angle=2 * PI - angle,
            start_angle=angle,
            color=BLUE,
            fill_opacity=0.5
        )
        self.add(lower_fill)
        lower_semicircle.add(lower_fill)
        self.remove(temp_circle)


        self.play(Create(upper_seperator_1))
        self.add(upper_seperator_2)
        arc_center = upper_semicircleArc.get_arc_center()
        self.play(
            Rotate(
                upper_seperator_2, 
                angle = angle,  
                about_point = arc_center
            ),
            run_time=1
        )
        self.add(lower_seperator)

        # Animate the group moving to the left
        # self.play(semi_circles.animate.shift(LEFT * 2))

        

        self.play(upper_semicircle.animate.shift(RIGHT * 2),
                lower_semicircle.animate.shift(LEFT * 2))
        self.wait(1)

        self.play(FadeOut(lower_semicircle))
        self.remove(lower_semicircle)

        self.play(upper_semicircle.animate.shift(LEFT * 2))

        # Set up camera for 3D viewing
        self.set_camera_orientation(phi=0 * DEGREES, theta=270 * DEGREES, zoom=1)
        self.move_camera(phi=75 * DEGREES, theta=30 * DEGREES, zoom=2.5, run_time=2)
        self.wait(0.5)

        # Create a copy of the semicircle to be deformed
        # deforming_semicircle = upper_semicircle.copy()
        # self.remove(upper_semicircle)  # Hide the original semicircle
        # self.add(deforming_semicircle)

        # Move the semicircle back to center first
        # self.play(upper_semicircle.animate.move_to(ORIGIN), run_time=1)

        
        self.play(Indicate(upper_semicircle.submobjects[0]), Indicate(upper_semicircle.submobjects[1]), run_time=3)


        # Create a copy of the semicircle to animate
        deforming_semicircle = upper_semicircle.copy()
        # self.add(deforming_semicircle)
        # self.remove(upper_semicircle)  # Hide the original

        # Define the animation
        def semicircle_to_cone(mob, alpha):
            # alpha goes from 0 to 1 during the animation
            
            # Get the original parts
            arc = mob.submobjects[2]  # The arc is the third submobject
            line1 = mob.submobjects[0]  # First separator line
            line2 = mob.submobjects[1]  # Second separator line
            
            # Original center and radius
            center = arc.get_arc_center()
            r = radius  # Radius stays 1
            
            # r = ((radius * (2*PI - angle * alpha)/(2*PI)))
            r = radius * (1 - alpha *(1 - angle/(2*PI)))  
            print(f"\nr:{r}  radius: {radius}  alpha: {alpha}")
            
            
            # Calculate height based on alpha
            height = np.sqrt(radius**2 - r**2)

            
            base_radius = r
            
            # Calculate the arc angle - as alpha increases, it approaches 2*PI (full circle)
            # arc_angle = angle + alpha * (2*PI - angle)
            arc_angle = angle * (1 - alpha) + 2*PI * alpha
            arc_start = 0  # Center the arc
            
            # Calculate endpoints for the arc and lines
            apex = center + height * OUT  # Top point of the cone
            
            # Calculate endpoints of the lines at the base of the cone
            # As alpha increases, the lines converge to form a full circle
            line1_base = center + base_radius * np.array([np.cos(arc_start), np.sin(arc_start), 0])
            line2_base = center + base_radius * np.array([np.cos(arc_start + arc_angle), np.sin(arc_start + arc_angle), 0])
            
            # Create the updated arc that is gradually becoming a circle
            new_arc = Arc(
                radius=base_radius,
                start_angle=arc_start,
                angle=arc_angle,
                color=BLUE,
                stroke_width=4  # Make arc more visible
            )
            new_arc.shift(center)  # Position at center
            
            # Update the lines to connect to the apex
            line1.set_points_as_corners([line1_base, apex])
            line2.set_points_as_corners([line2_base, apex])
            
            # Update the fill
            fill = mob.submobjects[3]  # The fill is the fourth submobject

            num_points = 20
            
            # Create a proper 3D surface for the cone mantle
            # We'll use a parametric surface to create the curved surface
            z_height = height * (1 - np.linspace(0, 1, num_points+1))  # Heights decrease from apex to base
            radii = np.linspace(0, base_radius, num_points+1)  # Radii increase from apex to base

            # Create 3D points for the mantle surface
            surface_points = []
            for i in range(num_points + 1):
                # Generate points in a circle at each height
                r = radii[i]
                z = z_height[i]
                # For each height, add points along the arc
                for j in range(num_points + 1):
                    t = j / num_points
                    angle_param = arc_start + t * arc_angle
                    x = center[0] + r * np.cos(angle_param)
                    y = center[1] + r * np.sin(angle_param)
                    surface_points.append(np.array([x, y, center[2] + z]))

            # Create polygons connecting the points to form the surface
            polys = []
            for i in range(num_points):
                for j in range(num_points):
                    # Get the four corners of each quad
                    i_j = i * (num_points + 1) + j
                    i1_j = (i + 1) * (num_points + 1) + j
                    i_j1 = i * (num_points + 1) + (j + 1)
                    i1_j1 = (i + 1) * (num_points + 1) + (j + 1)
                    
                    # Create two triangles for each quad
                    poly1 = Polygon(
                        surface_points[i_j], surface_points[i1_j], surface_points[i_j1],
                        color=BLUE, fill_opacity=0.5, stroke_width=0
                    )
                    poly2 = Polygon(
                        surface_points[i1_j], surface_points[i1_j1], surface_points[i_j1],
                        color=BLUE, fill_opacity=0.5, stroke_width=0
                    )
                    polys.append(poly1)
                    polys.append(poly2)

            # Create a VGroup with all the surface polygons
            new_fill = VGroup(*polys)

            # Replace the old fill
            mob.submobjects[3] = new_fill
            
            # Replace the old arc
            mob.submobjects[2] = new_arc
            


        # Create and play the animation
        self.remove(upper_semicircle)  # Hide the original semicircle
        self.play(
            UpdateFromAlphaFunc(deforming_semicircle, semicircle_to_cone),
            run_time=3
        )
        self.wait(1)

        # Create the final cone once the animation is complete
        cone = Cone(
            height= np.sqrt(radius**2 - (radius*(angle)/(2*PI))**2),
            base_radius=angle/(2*PI),
            direction=OUT,
            color=BLUE,
            fill_opacity=0.5,
            show_base=True
        )

        cone.move_to(deforming_semicircle.get_center())

        # Replace the deformed semicircle with the clean cone
        # self.remove(deforming_semicircle)
        self.add(cone)
        self.remove(deforming_semicircle)

        # Rotate camera to show the 3D cone
        self.move_camera(phi=90 * DEGREES, theta=0 * DEGREES, zoom=2.5, run_time=2)
        self.wait(0.5)

        arc = Arc(radius=1, start_angle=0, angle=2*PI)
        arc.set_color(BLUE)
        arc.set_stroke(width=2)
        seperator_1 = Line(start=arc.get_start(), end=arc.get_arc_center())
        seperator_2 = Line(start=arc.get_end(), end=arc.get_arc_center())
        seperator_1.set_color(WHITE)
        seperator_2.set_color(WHITE)
        seperator_1.set_stroke(width=1)
        seperator_2.set_stroke(width=1)

        fill = AnnularSector(
            inner_radius=0,
            outer_radius=1,
            angle=angle,
            start_angle=0,
            color=BLUE,
            fill_opacity=0.5
        )

        semi_circle = VGroup(arc, seperator_1, seperator_2, fill)
        semi_circle.rotate(PI / 2, axis=UP)  # Rotate to make it parallel to the camera

        semi_circle.shift(DOWN * 4)  # Move it up to match the cone's position

        self.play(semi_circle.animate.shift(UP * 3), cone.animate.shift(UP * 2), run_time=2)

        self.wait(1)

        # Animate the semi_circles angle going up to 270 degrees and back down to the original angle
        new_angle = 3 * PI / 2  # 270 degrees in radians

        def update_angle_semi_circle(mob, alpha):
            current_angle = angle + alpha * (new_angle - angle)
            mob.submobjects[0].become(Arc(radius=1, start_angle=0, angle=current_angle, color=BLUE))
            mob.submobjects[1].become(Line(start=mob.submobjects[0].get_start(), end=mob.submobjects[0].get_arc_center()))
            mob.submobjects[2].become(Line(start=mob.submobjects[0].get_end(), end=mob.submobjects[0].get_arc_center()))
            mob.submobjects[3].become(AnnularSector(
                inner_radius=0,
                outer_radius=1,
                angle=current_angle,
                start_angle=0,
                color=BLUE,
                fill_opacity=0.5
            ))
            mob.rotate(PI / 2, axis=UP)  # Ensure it keeps its rotation
            mob.shift(DOWN * 1)  # Ensure it keeps its position

        def update_angle_cone(mob, alpha):
            current_angle = angle + alpha * (new_angle - angle)
            # Recreate the cone with new dimensions
            new_cone = Cone(
                height=np.sqrt(radius**2 - (radius*(current_angle)/(2*PI))**2),
                base_radius=radius*(current_angle)/(2*PI),
                direction=OUT,
                color=BLUE,
                fill_opacity=0.5,
                show_base=True
            )
            new_cone.move_to(mob.get_center())
            mob.become(new_cone)
            # mob.shift(DOWN * 2)  # Ensure it keeps its position

        self.play(UpdateFromAlphaFunc(semi_circle, update_angle_semi_circle), UpdateFromAlphaFunc(cone, update_angle_cone), run_time=2)
        self.wait(0.5)

        def reverse_angle_semi_circle(mob, alpha):
            current_angle = new_angle - alpha * (new_angle - angle)
            mob.submobjects[0].become(Arc(radius=1, start_angle=0, angle=current_angle, color=BLUE))
            mob.submobjects[1].become(Line(start=mob.submobjects[0].get_start(), end=mob.submobjects[0].get_arc_center()))
            mob.submobjects[2].become(Line(start=mob.submobjects[0].get_end(), end=mob.submobjects[0].get_arc_center()))
            mob.submobjects[3].become(AnnularSector(
            inner_radius=0,
            outer_radius=1,
            angle=current_angle,
            start_angle=0,
            color=BLUE,
            fill_opacity=0.5
            ))
            mob.rotate(PI / 2, axis=UP)  # Ensure it keeps its rotation
            mob.shift(DOWN * 1)  # Ensure it keeps its position

        def reverse_angle_cone(mob, alpha):
            current_angle = new_angle - alpha * (new_angle - angle)
            # Recreate the cone with new dimensions
            new_cone = Cone(
                height=np.sqrt(radius**2 - (radius*(current_angle)/(2*PI))**2),
                base_radius=radius*(current_angle)/(2*PI),
                direction=OUT,
                color=BLUE,
                fill_opacity=0.5,
                show_base=True
            )
            new_cone.move_to(mob.get_center())
            mob.become(new_cone)
            # mob.shift(DOWN * 2)  # Ensure it keeps its position



        self.play(UpdateFromAlphaFunc(semi_circle, reverse_angle_semi_circle), UpdateFromAlphaFunc(cone, reverse_angle_cone), run_time=2)
        self.wait(0.5)




        self.wait(1)



class part2(ThreeDScene):
    def construct(self):
        
        self.set_camera_orientation(gamma=180 * DEGREES, phi=0 * DEGREES, theta=90 * DEGREES, zoom=2)
        cone = Cone(
            height= np.sqrt(radius**2 - (radius*(angle)/(2*PI))**2),
            base_radius=angle/(2*PI),
            direction=IN,
            color=BLUE,
            fill_opacity=0.2,
            show_base=True
            )
        self.add(cone)

        cone.move_to(ORIGIN)  # Center the cone in the scene
        cone.rotate(PI / 3, axis=RIGHT)  # Rotate to make it parallel to the camera
        cone.shift(RIGHT * 2)  # Move it up to match the cone's position

        self.play(cone.animate.shift(LEFT * 4.5), run_time=2)

        cone_formula = MathTex(
            r"V = \frac{\pi r^2 h}{3}"
        )
        height_line = Line(start=cone.get_center(), end=cone.get_top())
        self.play(Write(cone_formula))
        self.add_fixed_in_frame_mobjects(cone_formula)
        self.play(Indicate(cone_formula[0][-3]), run_time=2)  # Indicate only the "h"
        self.play(Create(height_line))
        self.play(Indicate(height_line), Indicate(cone_formula[0][-3]), run_time=2)
        
        self.wait(2)

        self.play(Indicate(cone_formula[0][-4]), Indicate(cone_formula[0][-5]), run_time=2)
        radius_line = Line(start=cone.get_center(), end=cone.get_center() + cone.width/2 * RIGHT)
        self.play(Create(radius_line))
        self.play(Indicate(radius_line), Indicate(cone_formula[0][-4]), Indicate(cone_formula[0][-5]), run_time=2)
        

        cone_formula.shift(UP * 2)

        radius_formula = MathTex(
            r"R(O) = \frac{O}{2\pi}" 
        )
        self.play(Write(radius_formula), run_time=2)
        self.add_fixed_in_frame_mobjects(radius_formula)

        self.play(Indicate(radius_formula[0][-4]), run_time=2)  # Indicate only the "R"
        self.play(radius_formula.animate.shift(UP * 2.5, RIGHT * 3.5), cone_formula.animate.shift(UP * 0.5))

        self.wait(5)

        seconfrius_formula = MathTex(
            r"O(v) = \frac{2\pi s(2\pi-v)}{2\pi}"
        )
        self.play(Write(seconfrius_formula), run_time=2)
        self.add_fixed_in_frame_mobjects(seconfrius_formula)

        self.wait(5)

        self.play(seconfrius_formula.animate.shift(UP * 1.5, RIGHT * 3.5))

        self.play(cone.animate.shift(UP * 1.5), height_line.animate.shift(UP * 1.5), radius_line.animate.shift(UP * 1.5), run_time=2)


        temp_circle = Circle(radius=1, color=BLUE, fill_opacity=0.5)
        self.play(Create(temp_circle))
        


        upper_semicircleArc = Arc(radius=1, start_angle=0, angle=angle)
        lower_semicircleArc = Arc(radius=1, start_angle=angle, angle=2*PI - angle)

        upper_semicircle = VGroup(upper_semicircleArc)
        lower_semicircle = VGroup(lower_semicircleArc)

        # Create a group with the 2 semicircles
        semi_circles = VGroup(upper_semicircleArc, lower_semicircleArc)
        # Set the color of the semicircles
        semi_circles.set_color(BLUE)


        # Add the group to the scene
        self.add(upper_semicircleArc)
        self.add(lower_semicircleArc)
        self.wait(1)
        # split the lines into 2 to allow for other angles then 180 degrees
        upper_seperator_1 = Line(start=upper_semicircleArc.get_start(), end=upper_semicircleArc.get_arc_center())
        upper_seperator_2 = Line(start=upper_semicircleArc.get_start(), end=upper_semicircleArc.get_arc_center())

        lower_seperator_1 = Line(start=lower_semicircleArc.get_start(), end=lower_semicircleArc.get_arc_center())
        lower_seperator_2 = Line(start=lower_semicircleArc.get_arc_center(), end=lower_semicircleArc.get_end())  
        
        upper_seperator = VGroup(upper_seperator_1, upper_seperator_2)
        lower_seperator = VGroup(lower_seperator_1, lower_seperator_2)

        upper_semicircle.add(upper_seperator_1, upper_seperator_2, upper_semicircleArc)
        lower_semicircle.add(lower_seperator_1, lower_seperator_2, lower_semicircleArc)

        # Fill the upper semicircle
        upper_fill = AnnularSector(
            inner_radius=0,
            outer_radius=1,
            angle=angle,
            start_angle=0,
            color=BLUE,
            fill_opacity=0.5
        )
        self.add(upper_fill)
        upper_semicircle.add(upper_fill)

        # Fill the lower semicircle
        lower_fill = AnnularSector(
            inner_radius=0,
            outer_radius=1,
            angle=2 * PI - angle,
            start_angle=angle,
            color=BLUE,
            fill_opacity=0.5
        )
        self.add(lower_fill)
        lower_semicircle.add(lower_fill)
        self.remove(temp_circle)

        self.play(Create(upper_seperator_1))
        self.add(upper_seperator_2)
        arc_center = upper_semicircleArc.get_arc_center()
        self.play(
            Rotate(
                upper_seperator_2, 
                angle = angle,  
                about_point = arc_center
            ),
            run_time=1
        )
        self.add(lower_seperator)
        circle = VGroup(semi_circles, upper_semicircle, lower_semicircle, upper_seperator_1, upper_seperator_2, lower_seperator_1, lower_seperator_2)

        self.play(circle.animate.shift(LEFT * 2.5))

        # Shrink the circle to half its size
        self.play(circle.animate.scale(0.5))




        self.wait(2)
        self.play(Indicate(upper_semicircleArc), Indicate(seconfrius_formula[0][-5]), run_time=3)

        self.play(Indicate(upper_seperator_1), Indicate(seconfrius_formula[0][-10]), run_time=3)

        

        


        




class OpeningManimExample(Scene):
    def construct(self):
        intro_words = Text("""
            The original motivation for manim was to
            better illustrate mathematical functions
            as transformations.
        """)
        intro_words.to_edge(UP)

        self.play(Write(intro_words))
        self.wait(2)

        # Linear transform
        grid = NumberPlane(x_range=(-10, 10), y_range=(-5, 5))
        matrix = [[1, 1], [0, 1]]
        linear_transform_words = VGroup(
            Text("This is what the matrix"),
            Matrix(matrix),
            Text("looks like")
        )
        linear_transform_words.arrange(RIGHT)
        linear_transform_words.to_edge(UP)

        self.play(
            Create(grid),
            FadeTransform(intro_words, linear_transform_words)
        )
        self.wait()
        self.play(grid.animate.apply_matrix(matrix), run_time=3)
        self.wait()

        # Complex map
        c_grid = ComplexPlane()
        moving_c_grid = c_grid.copy()
        moving_c_grid.prepare_for_nonlinear_transform()
        c_grid.set_stroke(BLUE_E, 1)
        c_grid.add_coordinates()
        complex_map_words = MathTex(r"""
            \text{Or thinking of the plane as } \mathbb{C},\\
            \text{this is the map } z \rightarrow z^2
        """)
        complex_map_words.to_corner(UR)

        self.play(
            FadeOut(grid),
            Create(c_grid, run_time=3),
            FadeIn(moving_c_grid),
            FadeTransform(linear_transform_words, complex_map_words),
        )
        self.wait()
        self.play(
            moving_c_grid.animate.apply_complex_function(lambda z: z**2),
            run_time=6,
        )
        self.wait(2)