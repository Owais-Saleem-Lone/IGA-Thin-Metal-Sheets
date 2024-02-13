import OCC
from OCC.Display.SimpleGui import init_display
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

def create_geometry(display):
    # Your code to create or load 3D geometry
    # For example, create a box
    box = BRepPrimAPI_MakeBox(10., 20., 30.).Shape()
    
    # Display the geometry
    display.DisplayShape(box, update=True)

def on_edge_selected():
    print("Edge selected")

def on_face_selected():
    print("Face selected")

if __name__ == '__main__':
    display, start_display, add_menu, add_function_to_menu = init_display()

    # Create a 'Selection' menu with custom functions
    add_menu('Selection')
    add_function_to_menu('Selection', on_edge_selected)
    add_function_to_menu('Selection', on_face_selected)

    # Create and display 3D geometry (the box)
    create_geometry(display)

    # Start the interactive display
    start_display()