import rhinoscriptsyntax as rs
import scriptcontext as sc
import Rhino

def get_trimmed_surface_details(trimmed_surface_id):
    # Get the trimming curves of the trimmed surface
    trim_curves = rs.DuplicateSurfaceBorder(trimmed_surface_id, type=2)

    if trim_curves:
        # Get the trimmed surface
        trimmed_surface = sc.doc.Objects.Find(trimmed_surface_id).Geometry

        # Print details of the trimming curves in parametric space
        for i, trim_curve_id in enumerate(trim_curves):
            trim_curve = rs.coercecurve(trim_curve_id, -1, True)
            if trim_curve:
                # Get the parameter space curve
                u_domain = trim_curve.Domain(0)
                v_domain = trim_curve.Domain(1)

                # Sample the parametric curve
                u_values = [u_domain.ParameterAt(t) for t in u_domain.GetBasisBlind(0)]
                v_values = [v_domain.ParameterAt(t) for t in v_domain.GetBasisBlind(0)]

                # Print parametric space details
                print(f"Parametric Space Details for Trim Curve {i + 1}:")
                print("U Values:", u_values)
                print("V Values:", v_values)
                print("\n")

                # Clean up temporary curve
                rs.DeleteObject(trim_curve_id)

def main():
    # Replace 'your_trimmed_surface_id' with the actual ID of your trimmed surface
    trimmed_surface_id = rs.GetObject("Select the trimmed surface", rs.filter.surface)

    if trimmed_surface_id:
        get_trimmed_surface_details(trimmed_surface_id)
    else:
        print("No trimmed surface selected.")

if __name__ == "__main__":
    main()
