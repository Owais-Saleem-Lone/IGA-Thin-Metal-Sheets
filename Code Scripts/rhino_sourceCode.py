import rhinoscriptsyntax as rs

def get_surface_control_points(surface_id):
    if rs.IsSurface(surface_id):
        control_points = rs.SurfaceEditPoints(surface_id)
        weights = rs.SurfaceWeights(surface_id)

        if control_points and weights:
            return control_points, weights
        else:
            print("Error: Unable to retrieve control points and weights.")
    else:
        print("Error: Input is not a valid surface.")

# Example usage:
# Replace 'your_surface_id' with the actual ID of your surface
surface_id = rs.GetObject("Select a surface", rs.filter.surface)
if surface_id:
    control_points, weights = get_surface_control_points(surface_id)

    if control_points and weights:
        print("Control Points:")
        for i, point in enumerate(control_points):
            print(f"Point {i+1}: {point}")

        print("\nWeights:")
        for i, weight in enumerate(weights):
            print(f"Weight {i+1}: {weight}")
