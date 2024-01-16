# Global Imports
import numpy as np
import FlowCal
from matplotlib.path import Path

def _gate_points(
    flowcal_object: FlowCal.io.FCSData,
    channels: list,
    gate_vertices: np.array
    ) -> np.array:
    """
    Args:
    flowcal_object: FCSData object to be gated
    channels: The channel axes of the plot to be gated
    gate_vertices: The vertices of the gate as a numpy array of shape (n, 2)
    
    Returns a numpy mask based on whether or not the unnown points are within the polygon defined by given vertices"""
    path = Path(gate_vertices)
    
    # Select only channels of interest (the two dimensions concerned with this particular gate)
    point_cloud = flowcal_object[:, channels]

    # determine if each point in unknown_array is within the polygon
    in_polygon = path.contains_points(point_cloud)
    return in_polygon

def gate_contour(
    flowcal_object: FlowCal.io.FCSData,
    channels: list,
    contour: list,
    ) -> np.array:
    """
    Args:
    flowcal_object: FCSData object to be gated\n
    channels: The channel axes of the plot to be gated\n
    contour: A list of gate vertices (which are numpy arrays of shape (n, 2))\n
    
    Returns a numpy mask based on whether or not the unnown points are within any polygons in the contour list"""
    
    contour_list = []
    for polygon in contour:
        contour_list.append(
            _gate_points(
                flowcal_object=flowcal_object,
                channels=channels,
                gate_vertices=polygon
            )
        )    
    return np.any(contour_list, axis=0)
    
def scale_polygon(
    vertices: np.array,
    scale_factor: float
    ) -> np.array:
    """
    Args:
    vertices: The vertices of the polygon in question as a numpy array of shape (n, 2)
    scale_factor: The factor by which to enlarge the polygon. For example: 1.2 will result in a 20% increase in the distance between the polygons euclidean centroid and each of its vertices. If you want to scale according to area instead, the best naiive approximation would be to multiply by 0.56419 (or the square root of 1/pi).
    """
    centroid = np.mean(vertices, axis=0)
    translated_vertices = vertices - centroid
    scaled_vertices = translated_vertices * scale_factor
    realigned_vertices = scaled_vertices + centroid

    return realigned_vertices
