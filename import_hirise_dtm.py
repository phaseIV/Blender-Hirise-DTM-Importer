# This file is a part of the HiRISE DTM Importer for Blender
#
# Copyright (C) 2017 Arizona Board of Regents on behalf of the Planetary Image
# Research Laboratory, Lunar and Planetary Laboratory at the University of
# Arizona.
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""A HiRISE DTM Importer for Blender"""

import bpy

import bpy.props
from bpy_extras.io_utils import ImportHelper
import bmesh

import numpy as np

import re
import collections


bl_info = {
    "name": "HiRISE DTM Importer",
    "author": "Nicholas Wolf (nicwolf@pirl.lpl.arizona.edu) / phaseIV",
    "version": (0, 2, 0),
    "blender": (2, 80, 0),
    "license": "GPL",
    "location": "File > Import > HiRISE DTM (.img)",
    "description": "Import a HiRISE DTM as a mesh.",
    "warning": "May consume a lot of memory",
    "category": "Import-Export",
    "wiki_url": "",  # TBD
    "tracker_url": "",  # TBD
    "link": "",  # TBD
    "support": "TESTING",
}

'''
if "bpy" in locals():
    import imp
    imp.reload(importer)
    imp.reload(terrainpanel)
'''


#------- ui/importer.py


class ImportHiRISETerrain(bpy.types.Operator, ImportHelper):
    """DTM Import Helper"""
    bl_idname = "import.pds_dtm"
    bl_label = "Import HiRISE Terrain Model"
    bl_options = {'UNDO'}
    filename_ext = ".img"
    filter_glob: bpy.props.StringProperty(
        options={'HIDDEN'},
        default="*.img"
    )

    # Allow the user to specify a resolution factor for loading the
    # terrain data at. This is useful because it allows the user to stage
    # a scene with a low resolution terrain map, apply textures, modifiers,
    # etc. and then increase the resolution to prepare for final rendering.
    #
    # Displaying this value as a percentage (0, 100] is an intuitive way
    # for users to grasp what this value does. The DTM importer, however,
    # wants to recieve a value between (0, 1]. This is obviously a
    # straightforward conversion:
    #
    #     f(x) = x / 100
    #
    # But this conversion should happen here, in the terrain panel, rather
    # than in the DTM importing utility itself. We can't pass get/set
    # functions to the property itself because they result in a recursion
    # error. Instead, we use another, hidden, property to store the scaled
    # resolution.
    dtm_resolution: bpy.props.FloatProperty(
        subtype="PERCENTAGE",
        description=(
            "Percentage scale for terrain model resolution. 100\% loads the "
            "model at full resolution (i.e. one vertex for each post in the "
            "original terrain model) and is *MEMORY INTENSIVE*. Downsampling "
            "uses Nearest Neighbors. You will be able to increase the "
            "resolution of your mesh later, and still maintain all textures, "
            "transformations, modifiers, etc., so best practice is to start "
            "small. The downsampling algorithm may need to alter the "
            "resolution you specify here to ensure it results in a whole "
            "number of vertices. If it needs to alter the value you specify, "
            "you are guaranteed that it will shrink it (i.e. decrease the "
            "DTM resolution."
        ),
        name="Terrain Model Resolution",
        min=1.0, max=100.0, default=10.0
    )
    scaled_dtm_resolution: bpy.props.FloatProperty(
        options={'HIDDEN'},
        name="Scaled Terrain Model Resolution",
        get=(lambda self: self.dtm_resolution / 100)
    )

    # HiRISE DTMs are huge, but it can be nice to load them in at scale. Here,
    # we present the user with the option of setting up the Blender viewport
    # to avoid a couple of common pitfalls encountered when working with such
    # a large mesh.
    #
    # 1. The Blender viewport has a default clipping distance of 1km. HiRISE
    #    DTMs are often many kilometers in each direction. If this setting is
    #    not changed, an unsuspecting user may only see part (or even nothing
    #    at all) of the terrain. This option (true, by default) instructs
    #    Blender to change the clipping distance to something appropriate for
    #    the DTM, and scales the grid floor to have gridlines 1km apart,
    #    instead of 1m apart.
    should_setup_viewport: bpy.props.BoolProperty(
        description=(
            "Set up the Blender screen to try and avoid clipping the DTM "
            "and to make the grid floor larger. *WARNING* This will change "
            "clipping distances and the Blender grid floor, and will fit the "
            "DTM in the viewport."
        ),
        name="Setup Blender Scene", default=False
    )
    # 2. Blender's default units are dimensionless. This option instructs
    #    Blender to change its unit's dimension to meters.
    should_setup_units: bpy.props.BoolProperty(
        description=(
            "Set the Blender scene to use meters as its unit."
        ),
        name="Set Blender Units to Meters", default=False
    )

    def execute(self, context):
        """Runs when the "Import HiRISE Terrain Model" button is pressed"""
        filepath = bpy.path.ensure_ext(self.filepath, self.filename_ext)
        # Create a BTerrain from the DTM
        dtm = DTM(filepath, self.scaled_dtm_resolution)
        BTerrain.new(dtm)

        # Set up the Blender UI
        if self.should_setup_units:
            self._setup_units(context)
        if self.should_setup_viewport:
            self._setup_viewport(context)

        return {"FINISHED"}

    def _setup_units(self, context):
        """Sets up the Blender scene for viewing the DTM"""
        scene = bpy.context.scene

        # Set correct units
        scene.unit_settings.system = 'METRIC'
        scene.unit_settings.scale_length = 1.0

        return {'FINISHED'}

    def _setup_viewport(self, context):
        """Sets up the Blender screen to make viewing the DTM easier"""
        screen = bpy.context.screen
        '''
        # Fetch the 3D_VIEW Area
        for area in screen.areas:
            if area.type == 'VIEW_3D':
                space = area.spaces[0]
                # Adjust 3D View Properties
                # TODO: Can these be populated more intelligently?
                space.clip_end = 100000
                space.grid_scale = 1000
                space.grid_lines = 50
        '''
        # Fly to a nice view of the DTM
        self._view_dtm(context)

        return {'FINISHED'}

    def _view_dtm(self, context):
        """Sets up the Blender screen to make viewing the DTM easier"""
        screen = bpy.context.screen

        # Fetch the 3D_VIEW Area
        for area in screen.areas:
            if area.type == 'VIEW_3D':
                # Move the camera around in the viewport. This requires
                # a context override.
                for region in area.regions:
                    if region.type == 'WINDOW':
                        override = {
                            'area': area,
                            'region': region,
                            'edit_object': bpy.context.edit_object
                        }
                        # Center View on DTM (SHORTCUT: '.')
                        bpy.ops.view3d.view_selected(override)
                        # Move to 'TOP' viewport (SHORTCUT: NUMPAD7)
                        #bpy.ops.view3d.viewnumpad(override, type='TOP')

        return {'FINISHED'}


#------- ui/terrainpanel.py


class TerrainPanel(bpy.types.Panel):
    """Creates a Panel in the Object properites window for terrain objects"""
    bl_label = "Terrain Model"
    bl_idname = "OBJECT_PT_terrain"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "object"

    # Allow the user to specify a new resolution factor for reloading the
    # terrain data at. This is useful because it allows the user to stage
    # a scene with a low resolution terrain map, apply textures, modifiers,
    # etc. and then increase the resolution to prepare for final rendering.
    #
    # Displaying this value as a percentage (0, 100] is an intuitive way
    # for users to grasp what this value does. The DTM importer, however,
    # wants to recieve a value between (0, 1]. This is obviously a
    # straightforward conversion:
    #
    #     f(x) = x / 100
    #
    # But this conversion should happen here, in the terrain panel, rather
    # than in the DTM importing utility itself. We can't pass get/set
    # functions to the property itself because they result in a recursion
    # error. Instead, we use another, hidden, property to store the scaled
    # resolution.
    bpy.types.Object.dtm_resolution = bpy.props.FloatProperty(
        subtype="PERCENTAGE",
        name="New Resolution",
        description=(
            "Percentage scale for terrain model resolution. 100\% loads the "
            "model at full resolution (i.e. one vertex for each post in the "
            "original terrain model) and is *MEMORY INTENSIVE*. Downsampling "
            "uses Nearest Neighbors. The downsampling algorithm may need to "
            "alter the resolution you specify here to ensure it results in a "
            "whole number of vertices. If it needs to alter the value you "
            "specify, you are guaranteed that it will shrink it (i.e. "
            "decrease the DTM resolution."
        ),
        min=1.0, max=100.0, default=10.0
    )
    bpy.types.Object.scaled_dtm_resolution = bpy.props.FloatProperty(
        options={'HIDDEN'},
        name="Scaled Terrain Model Resolution",
        get=(lambda self: self.dtm_resolution / 100.0)
    )

    @classmethod
    def poll(cls, context):
        return context.object.get("IS_TERRAIN", False)
    
    def draw(self, context):
        obj = context.active_object
        layout = self.layout

        # User Controls
        layout.prop(obj, 'dtm_resolution')
        layout.operator("terrain.reload")

        # Metadata
        self.draw_metadata_panel(context)

    def draw_metadata_panel(self, context):
        """Display some metadata about the DTM"""
        obj = context.active_object
        layout = self.layout

        metadata_panel = layout.box()

        dtm_resolution = metadata_panel.row()
        dtm_resolution.label(text='Current Resolution: ')
        dtm_resolution.label(text='{:9,.2%}'.format(
            obj['DTM_RESOLUTION']
        ))

        mesh_scale = metadata_panel.row()
        mesh_scale.label(text='Current Scale: ')
        mesh_scale.label(text='{:9,.2f} m/post'.format(
            obj['MESH_SCALE']
        ))

        dtm_scale = metadata_panel.row()
        dtm_scale.label(text='Original Scale: ')
        dtm_scale.label(text='{:9,.2f} m/post'.format(
            obj['MAP_SCALE']
        ))

        dtm_min_lat = metadata_panel.row()
        dtm_min_lat.label(text='Minimum Latitude: ')
        dtm_min_lat.label(text='{:9,.3f} lat'.format(
            obj['MINIMUM_LATITUDE']
        ))

        dtm_max_lat = metadata_panel.row()
        dtm_max_lat.label(text='Maximum Latitude: ')
        dtm_max_lat.label(text='{:9,.3f} lat'.format(
            obj['MAXIMUM_LATITUDE']
        ))

        dtm_east_lon = metadata_panel.row()
        dtm_east_lon.label(text='Easternmost Longitude: ')
        dtm_east_lon.label(text='{:9,.3f} lon'.format(
            obj['EASTERNMOST_LONGITUDE']
        ))

        dtm_west_lon = metadata_panel.row()
        dtm_west_lon.label(text='Westernmost Longitude: ')
        dtm_west_lon.label(text='{:9,.3f} lon'.format(
            obj['WESTERNMOST_LONGITUDE']
        ))

        return {'FINISHED'}


class ReloadTerrain(bpy.types.Operator):
    """Button for reloading the terrain mesh at a new resolution."""
    bl_idname = "terrain.reload"
    bl_label = "Reload Terrain"

    def execute(self, context):
        # Reload the terrain
        obj = context.object
        path = obj['PATH']

        scaled_dtm_resolution = obj.scaled_dtm_resolution

        # Reload BTerrain with new DTM
        print("Reload BTerrain with new DTM**")
        dtm = DTM(path, scaled_dtm_resolution)
        BTerrain.reload(obj, dtm)

        return {"FINISHED"}


#------- pvl/__init__.py


def PVLload(path):
    """Returns a dict-like representation of a PVL label"""
    return LabelParser.load(path)


#------- terrain.py


class BTerrain:
    """
    Functions for creating Blender meshes from DTM objects

    This class contains functions that convert DTM objects to Blender meshes.
    Its main responsiblity is to triangulate a mesh from the elevation data in
    the DTM. Additionally, it attaches some metadata to the object and creates
    a UV map for it so that companion ortho-images drape properly.

    This class provides two public methods: `new()` and `reload()`.

    `new()` creates a new object[1] and attaches a new mesh to it.

    `reload()` replaces the mesh that is attached to an already existing
    object. This allows us to retain the location and orientation of the parent
    object's coordinate system but to reload the terrain at a different
    resolution.

    Notes
    ----------
    [1] If you're unfamiliar with Blender, one thing that will help you in
        reading this code is knowing the difference between 'meshes' and
        'objects'. A mesh is just a collection of vertices, edges and
        faces. An object may have a mesh as a child data object and
        contains additional information, e.g. the location and orientation
        of the coordinate system its child-meshes are reckoned in terms of.

    """

    @staticmethod
    def new(dtm, name='Terrain'):
        """
        Loads a new terrain

        Parameters
        ----------
        dtm : DTM
        name : str, optional
            The name that will be assigned to the new object, defaults
            to 'Terrain' (and, if an object named 'Terrain' already
            exists, Blender will automatically extend the name of the
            new object to something like 'Terrain.001')

        Returns
        ----------
        obj : bpy_types.Object

        """
        bpy.ops.object.add(type="MESH")
        obj = bpy.context.object
        obj.name = name

        # Fill the object data with a Terrain mesh
        obj.data = BTerrain._mesh_from_dtm(dtm)

        # Add some meta-information to the object
        metadata = BTerrain._create_metadata(dtm)
        BTerrain._setobjattrs(obj, **metadata)

        # Center the mesh to its origin and create a UV map for draping
        # ortho images.
        BTerrain._center(obj)

        return obj

    @staticmethod
    def reload(obj, dtm):
        """
        Replaces an exisiting object's terrain mesh

        This replaces an object's mesh with a new mesh, transferring old
        materials over to the new mesh. This is useful for reloading DTMs
        at different resolutions but maintaining textures/location/rotation.

        Parameters
        -----------
        obj : bpy_types.Object
            An already existing Blender object
        dtm : DTM

        Returns
        ----------
        obj : bpy_types.Object

        """
        old_mesh = obj.data
        new_mesh = BTerrain._mesh_from_dtm(dtm)

        # Copy any old materials to the new mesh
        for mat in old_mesh.materials:
            new_mesh.materials.append(mat.copy())

        # Swap out the old mesh for the new one
        obj.data = new_mesh

        # Update out-dated meta-information
        metadata = BTerrain._create_metadata(dtm)
        BTerrain._setobjattrs(obj, **metadata)

        # Center the mesh to its origin and create a UV map for draping
        # ortho images.
        BTerrain._center(obj)

        return obj

    @staticmethod
    def _mesh_from_dtm(dtm, name='Terrain'):
        """
        Creates a Blender *mesh* from a DTM

        Parameters
        ----------
        dtm : DTM
        name : str, optional
            The name that will be assigned to the new mesh, defaults
            to 'Terrain' (and, if an object named 'Terrain' already
            exists, Blender will automatically extend the name of the
            new object to something like 'Terrain.001')

        Returns
        ----------
        mesh : bpy_types.Mesh

        Notes
        ----------
        * We are switching coordinate systems from the NumPy to Blender.

              Numpy:            Blender:
               + ----> (0, j)    ^ (0, y)
               |                 |
               |                 |
               v (i, 0)          + ----> (x, 0)

        """
        # Create an empty mesh
        mesh = bpy.data.meshes.new(name)

        # Get the xy-coordinates from the DTM, see docstring notes
        y, x = np.indices(dtm.data.shape).astype('float64')
        x *= dtm.mesh_scale
        y *= -1 * dtm.mesh_scale

        # Create an array of 3D vertices
        vertices = np.dstack([x, y, dtm.data]).reshape((-1, 3))

        # Drop vertices with NaN values (used in the DTM to represent
        # areas with no data)
        vertices = vertices[~np.isnan(vertices).any(axis=1)]

        
        # Calculate the faces of the mesh
        triangulation = Triangulate(dtm.data)
        faces = triangulation.face_list()

        # Fill the mesh
        mesh.from_pydata(vertices, [], faces)
        mesh.update()
        
        '''
        # Smooth the mesh
        for f in mesh.polygons:
            f.use_smooth = True
        '''
        
        # Create a new UV layer
        mesh.uv_layers.new(name="HiRISE Generated UV Map")

        # We'll use a bmesh to populate the UV map with values        
        bm = bmesh.new()
        bm.from_mesh(mesh)
        bm.faces.ensure_lookup_table()
        uv_layer = bm.loops.layers.uv[0]

        # Iterate over each face in the bmesh
        num_faces = len(bm.faces)
        w = dtm.data.shape[1]
        h = dtm.data.shape[0]
        for face_index in range(num_faces):
            # Iterate over each loop in the face
            for loop in bm.faces[face_index].loops:
                # Get this loop's vertex coordinates
                vert_coords = loop.vert.co.xy
                # And calculate it's uv coordinate. We do this by dividing the
                # vertice's x and y coordinates by:
                #
                #     d + 1, dimensions of DTM (in "posts")
                #     mesh_scale, meters/DTM "post"
                #
                # This has the effect of mapping the vertex to its
                # corresponding "post" index in the DTM, and then mapping
                # that value to the range [0, 1).
                u = vert_coords.x / ((w + 1) * dtm.mesh_scale)
                v = 1 + vert_coords.y / ((h + 1) * dtm.mesh_scale)
                loop[uv_layer].uv = (u, v)

        bm.to_mesh(mesh)
        return mesh

    @staticmethod
    def _center(obj):
        """Move object geometry to object origin"""
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.origin_set(center='BOUNDS')

    @staticmethod
    def _setobjattrs(obj, **attrs):
        for key, value in attrs.items():
            obj[key] = value

    @staticmethod
    def _create_metadata(dtm):
        """Returns a dict containing meta-information about a DTM"""
        return {
            'PATH': dtm.path,
            'MESH_SCALE': dtm.mesh_scale,
            'DTM_RESOLUTION': dtm.terrain_resolution,
            'BIN_SIZE': dtm.bin_size,
            'MAP_SIZE': dtm.map_size,
            'MAP_SCALE': dtm.map_scale * dtm.unit_scale,
            'UNIT_SCALE': dtm.unit_scale,
            'MINIMUM_LATITUDE': dtm.min_lat,
            'MAXIMUM_LATITUDE': dtm.max_lat,
            'EASTERNMOST_LONGITUDE': dtm.eastern_lon,
            'WESTERNMOST_LONGITUDE': dtm.western_lon,
            'IS_TERRAIN': True,
            'HAS_UV_MAP': True
        }


#------- dtm.py


class DTM:
    """
    HiRISE Digital Terrain Model

    This class imports a HiRISE DTM from a Planetary Data Systems (PDS)
    compliant .IMG file.

    Parameters
    ----------
    path : str
    terrain_resolution : float, optional
        Controls the resolution the DTM is read at. This should be a float
        in the range [0.01, 1.0] (and will be constrained to this range). A
        value of 1.0 will result in the DTM being read at full resolution. A
        value of 0.01 will result in the DTM being read at 1/100th resolution.
        Default is 1.0 (no downsampling).

    Todo
    ----
    * Use GDAL for importing the DTM if it is installed for this Python
      environment. If/when I have the time to do this, it probably
      warrants breaking out separate importer classes. The benefits of
      doing this are pretty substantial, though:

        + More reliable (doesn't rely on my PVL parser for finding the
          valid values in the DTM, for locating the starting position of
          the elevation data in the .IMG file)

        + Other, better, downsampling algorithms are already built in.

        + Would make this much better at general PDS DTM importing,
          currently some of the import code is specific to HiRISE DTMs.

    """

    # Special constants in our data:
    #     NULL : No data at this point.
    #     LRS  : Low Representation Saturation
    #     LIS  : Low Instrument Saturation
    #     HRS  : High Representation Saturation
    #     HIS  : High Insturment Saturation
    SPECIAL_VALUES = {
        "NULL": np.frombuffer(b'\xFF\x7F\xFF\xFB', dtype='>f4')[0],
        "LRS": np.frombuffer(b'\xFF\x7F\xFF\xFC', dtype='>f4')[0],
        "LIS": np.frombuffer(b'\xFF\x7F\xFF\xFD', dtype='>f4')[0],
        "HRS": np.frombuffer(b'\xFF\x7F\xFF\xFE', dtype='>f4')[0],
        "HIS": np.frombuffer(b'\xFF\x7F\xFF\xFF', dtype='>f4')[0]
    }

    def __init__(self, path, terrain_resolution=1.0):
        self.path = path
        self.terrain_resolution = terrain_resolution
        self.label = self._read_label()
        self.data = self._read_data()

    def _read_label(self):
        """Returns a dict-like representation of a PVL label"""
        return PVLload(self.path)

    def _read_data(self):
        """
        Reads elevation data from a PDS .IMG file.

        Notes
        -----
        * Uses nearest-neighbor to downsample data.

        Todo
        ----
        * Add other downsampling algorithms.

        """
        h, w = self.image_resolution
        max_samples = int(w - w % self.bin_size)

        data = np.zeros(self.shape)
        with open(self.path, 'rb') as f:
            # Seek to the first byte of data
            start_byte = self._get_data_start()
            f.seek(start_byte)
            # Iterate over each row of the data
            for r in range(data.shape[0]):
                # Each iteration, seek to the right location before
                # reading a row. We determine this location as the
                # first byte of data PLUS a offset which we calculate as the
                # product of:
                #
                #     4, the number of bytes in a single record
                #     r, the current row index
                #     w, the number of records in a row of the DTM
                #     bin_size, the number of records in a bin
                #
                # This is where we account for skipping over rows.
                offset = int(4 * r * w * self.bin_size)
                f.seek(start_byte + offset)
                # Read a row
                row = np.fromfile(f, dtype=np.float32, count=max_samples)
                # This is where we account for skipping over columns.
                data[r] = row[::self.bin_size]

        data = self._process_invalid_data(data)
        return data

    def _get_data_start(self):
        """Gets the start position of the DTM data block"""
        label_length = self.label['RECORD_BYTES']
        num_labels = self.label.get('LABEL_RECORDS', 1)
        return int(label_length * num_labels)

    def _process_invalid_data(self, data):
        """Sets any 'NULL' elevation values to np.NaN"""
        invalid_data_mask = (data <= self.SPECIAL_VALUES['NULL'])
        data[invalid_data_mask] = np.NaN
        return data

    @property
    def map_size(self):
        """Geographic size of the bounding box around the DTM"""
        scale = self.map_scale * self.unit_scale
        w = self.image_resolution[0] * scale
        h = self.image_resolution[1] * scale
        return (w, h)

    @property
    def mesh_scale(self):
        """Geographic spacing between mesh vertices"""
        return self.bin_size * self.map_scale * self.unit_scale

    @property
    def map_info(self):
        """Map Projection metadata"""
        return self.label['IMAGE_MAP_PROJECTION']

    @property
    def map_scale(self):
        """Geographic spacing between DTM posts"""
        map_scale = self.map_info.get('MAP_SCALE', None)
        return getattr(map_scale, 'value', 1.0)

    @property
    def min_lat(self):
        """MINIMUM_LATITUDE of DTM"""
        min_lat = self.map_info.get('MINIMUM_LATITUDE', None)
        return getattr(min_lat, 'value', 1.0)

    @property
    def max_lat(self):
        """MAXIMUM_LATITUDE of DTM"""
        max_lat = self.map_info.get('MAXIMUM_LATITUDE', None)
        return getattr(max_lat, 'value', 1.0)

    @property
    def eastern_lon(self):
        """EASTERNMOST_LONGITUDE of DTM"""
        eastern_lon = self.map_info.get('EASTERNMOST_LONGITUDE', None)
        return getattr(eastern_lon, 'value', 1.0)

    @property
    def western_lon(self):
        """WESTERNMOST_LONGITUDE of DTM"""
        western_lon = self.map_info.get('WESTERNMOST_LONGITUDE', None)
        return getattr(western_lon, 'value', 1.0)

    @property
    def map_units(self):
        """Geographic unit for spacing between DTM posts"""
        map_scale = self.map_info.get('MAP_SCALE', None)
        return getattr(map_scale, 'units', None)

    @property
    def unit_scale(self):
        """
        The function that creates a Blender mesh from this object will assume
        that the height values passed into it are in meters --- this
        property is a multiplier to convert DTM-units to meters.
        """
        scaling_factors = {
            'KM/PIXEL': 1000,
            'METERS/PIXEL': 1
        }
        return scaling_factors.get(self.map_units, 1.0)

    @property
    def terrain_resolution(self):
        """Vertex spacing, meters"""
        return self._terrain_resolution

    @terrain_resolution.setter
    def terrain_resolution(self, t):
        self._terrain_resolution = np.clip(t, 0.01, 1.0)

    @property
    def bin_size(self):
        """The width of the (square) downsampling bin"""
        return int(np.ceil(1 / self.terrain_resolution))

    @property
    def image_stats(self):
        """Image statistics from the original DTM label"""
        return self.label['IMAGE']

    @property
    def image_resolution(self):
        """(Line, Sample) resolution of the original DTM"""
        return (self.image_stats['LINES'], self.image_stats['LINE_SAMPLES'])

    @property
    def size(self):
        """Number of posts in our reduced DTM"""
        return self.shape[0] * self.shape[1]

    @property
    def shape(self):
        """Shape of our reduced DTM"""
        num_rows = self.image_resolution[0] // self.bin_size
        num_cols = self.image_resolution[1] // self.bin_size
        return (num_rows, num_cols)


#-------


class Triangulate:
    """
    A triangulation algorithm for creating a mesh from a DTM raster.

    I have been re-writing parts of the Blender HiRISE DTM importer in an
    effort to cull its dependencies on external packages. Originally, the
    add-on relied on SciPy's Delaunay triangulation (really a wrapper for
    Qhull's Delaunay triangulation) to triangulate a mesh from a HiRISE DTM.

    This re-write is much better suited to the problem domain. The SciPy
    Delaunay triangulation creates a mesh from any arbitrary point cloud and,
    while robust, doesn't care about the fact that our HiRISE DTMs are
    regularly gridded rasters. This triangulation algorithm is less robust
    but much faster. Credit is due to Tim Spriggs for his work on the previous
    Blender HiRISE DTM importer --- this triangulation algorithm largely
    models the one in his add-on with a few changes (namely interfacing
    with NumPy's API).

    Overview
    ----------
    Suppose we have a DTM:

    .. code::

                -  -  -  -  -  -  -  -  X  X  -  -  -  -  -
                -  -  -  -  -  -  X  X  X  X  X  -  -  -  -
                -  -  -  -  X  X  X  X  X  X  X  X  -  -  -
                -  -  X  X  X  X  X  X  X  X  X  X  X  -  -
                X  X  X  X  X  X  X  X  X  X  X  X  X  X  -
                -  X  X  X  X  X  X  X  X  X  X  X  X  X  X
                -  -  X  X  X  X  X  X  X  X  X  X  X  -  -
                -  -  -  X  X  X  X  X  X  X  X  -  -  -  -
                -  -  -  -  X  X  X  X  X  -  -  -  -  -  -
                -  -  -  -  -  X  X  -  -  -  -  -  -  -  -

    where 'X' represents valid values and '-' represents invalid values.
    Valid values should become vertices in the resulting mesh, invalid
    values should be ignored.

    Our end goal is to supply Blender with:

        1. an (n x 3) list of vertices

        2. an (m x 3) list of faces.

    A vertex is a 3-tuple that we get from the DTM raster array. The
    z-coordinate is whatever elevation value is in the DTM and the xy-
    coordinates are the image indices multiplied by the resolution of the
    DTM (e.g. if the DTM is at 5m/px, the first vertex is at (0m, 0m,
    z_00) and the vertex to the right of it is at (5m, 0m, z_01)).

    A face is a 3-tuple (because we're using triangles) where each element
    is an index of a vertex in the vertices list. Computing the faces is
    tricky because we want to leverage the orthogonal structure of the DTM
    raster for computational efficiency but we also need to reference
    vertex indices in our faces, which don't observe any regular
    structure.

    We take two rows at a time from the DTM raster and track the *raster
    row* indices as well as well as the *vertex* indices. Raster row
    indices are the distance of a pixel in the raster from the left-most
    (valid *or* invalid) pixel of the row. The first vertex is index 0 and
    corresponds to the upperleft-most valid pixel in the DTM raster.
    Vertex indices increase to the right and then down.

    For example, the first two rows:

    .. code::

                -  -  -  -  -  -  -  -  X  X  -  -  -  -  -
                -  -  -  -  -  -  X  X  X  X  X  -  -  -  -

    in vertex indices:

    .. code::

                -  -  -  -  -  -  -  -  0  1  -  -  -  -  -
                -  -  -  -  -  -  2  3  4  5  6  -  -  -  -

    and in raster row indices:

    .. code::

                -  -  -  -  -  -  -  -  9 10  -  -  -  -  -
                -  -  -  -  -  -  7  8  9 10 11  -  -  -  -

    To simplify, we will only add valid square regions to our mesh. So,
    for these first two rows the only region that will be added to our
    mesh is the quadrilateral formed by vertices 0, 1, 4 and 5. We
    further divide this area into 2 triangles and add the vertices to the
    face list in CCW order (i.e. t1: (4, 1, 0), t2: (4, 5, 1)).

    After the triangulation between two rows is completed, the bottom
    row is cached as the top row and the next row in the DTM raster is
    read as the new bottom row. This process continues until the entire
    raster has been triangulated.

    Todo
    ---------
    * It should be pretty trivial to add support for triangular
      regions (i.e. in the example above, also adding the triangles
      formed by (3, 4, 0) and (5, 6, 1)).

    """
    def __init__(self, array):
        self.array = array
        self.faces = self._triangulate()

    def _triangulate(self):
        """Triangulate a mesh from a topography array."""
        # Allocate memory for the triangles array
        max_tris = (self.array.shape[0] - 1) * (self.array.shape[1] - 1) * 2
        tris = np.zeros((max_tris, 3), dtype=int)
        ntri = 0

        # We initialize a vertex counter at 0
        prev_vtx_start = 0
        # We don't care about the values in the array, just whether or not
        # they are valid.
        prev = ~np.isnan(self.array[0])
        # We can sum this boolean array to count the number of valid entries
        prev_num_valid = prev.sum()
        # TODO: Probably a more clear (and faster) function than argmax for
        # getting the first Truth-y value in a 1d array.
        prev_img_start = np.argmax(prev)

        # Start quadrangulation
        for i in range(1, self.array.shape[0]):
            # Fetch this row, get our bearings in image *and* vertex space
            curr = ~np.isnan(self.array[i])
            curr_vtx_start = prev_vtx_start + prev_num_valid
            curr_img_start = np.argmax(curr)
            curr_num_valid = curr.sum()
            # Find the overlap between this row and the previous one
            overlap = np.logical_and(prev, curr)
            num_tris = overlap.sum() - 1
            overlap_start = np.argmax(overlap)
            # Store triangles
            for j in range(num_tris):
                curr_pad = overlap_start - curr_img_start + j
                prev_pad = overlap_start - prev_img_start + j
                tris[ntri + 0] = [
                    curr_vtx_start + curr_pad,
                    prev_vtx_start + prev_pad + 1,
                    prev_vtx_start + prev_pad
                ]
                tris[ntri + 1] = [
                    curr_vtx_start + curr_pad,
                    curr_vtx_start + curr_pad + 1,
                    prev_vtx_start + prev_pad + 1
                ]
                ntri += 2
            # Cache current row as previous row
            prev = curr
            prev_vtx_start = curr_vtx_start
            prev_img_start = curr_img_start
            prev_num_valid = curr_num_valid

        return tris[:ntri]

    def face_list(self):
        return list(self.faces)


#------- pvl/patterns.py


# End of PVL File
END = re.compile(
    r'\s* \bEND\b \s*', (re.VERBOSE + re.IGNORECASE)
)

# PVL Comment
COMMENT = re.compile(
    r'/\* .*? \*/', (re.DOTALL + re.VERBOSE)
)

# PVL Statement
STATEMENT = re.compile(
    r"""
    \s* (?P<key>\w+) # Match a PVL key
    \s* = \s* # Who knows how many spaces we encounter
    (?P<val> # Match a PVL value
        ([+-]?\d+\.?\d*) # We could match a number
        | (['"]?((\w+ \s*?)+)['"]?) # Or a string
    )
    (\s* <(?P<units>.*?) >)? # The value may have an associated unit
    """, re.VERBOSE
)

# Integer Number
INTEGER = re.compile(
    r"""
    [+-]?(?<!\.)\b[0-9]+\b(?!\.[0-9])
    """, re.VERBOSE
)

# Floating Point Number
FLOATING = re.compile(
    r"""
    [+-]?\b[0-9]*\.?[0-9]+
    """, re.VERBOSE
)


#------- pvl/label.py


class Label(dict):
    """A dict-like representation of a PVL label"""
    def __init__(self, *args, **kwargs):
        super(Label, self).__init__(*args, **kwargs)


#------- pvl/parse.py


Quantity = collections.namedtuple('Quantity', ['value', 'units'])


class PVLParseError(Exception):
    """Error parsing PVL file"""
    def __init__(self, message):
        super(PVLParseError, self).__init__(message)


class LabelParser:
    """A PVL Parser"""
    @staticmethod
    def load(path):
        """
        Load a dict-like representation of a PVL label header

        Parameters
        ----------
        path : str
            Path to a file containing a PVL header

        Returns
        ----------
        label : pvl.Label

        """
        raw = LabelParser._read(path)
        return Label(**LabelParser._parse(raw))

    @staticmethod
    def _read(path):
        """
        Get the PVL header from a file as a string

        Parameters
        ----------
        path : str
            Path to a file containing a PVL header

        Returns
        ----------
        raw : str

        Notes
        ---------
        * This function assumes that the file begins with a PVL header
          and it will read lines from the file until it encounters
          a PVL end statement.

        To-Do
        ---------
        * This could be more robust. What happens if there is no label
          in the file?

        """
        with open(path, 'rb') as f:
            raw = ''
            while True:
                try:
                    line = f.readline().decode()
                    raw += line
                    if re.match(END, line):
                        break
                except UnicodeDecodeError:
                    raise PVLParseError("Error parsing PVL label from "
                                        "file: {}".format(path))
        return raw

    @staticmethod
    def _remove_comments(raw):
        return re.sub(COMMENT, '', raw)

    @staticmethod
    def _parse(raw):
        raw = LabelParser._remove_comments(raw)
        label_iter = re.finditer(STATEMENT, raw)
        return LabelParser._parse_iter(label_iter)

    @staticmethod
    def _parse_iter(label_iter):
        """Recursively parse a PVL label iterator"""
        obj = {}
        while True:
            try:
                # Try to fetch the next match from the iter
                match = next(label_iter)
                val = match.group('val')
                key = match.group('key')
                # Handle nested object groups
                if key == 'OBJECT':
                    obj.update({
                        val: LabelParser._parse_iter(label_iter)
                    })
                elif key == 'END_OBJECT':
                    return obj
                # Add key/value pair to dict
                else:
                    # Should this value be a numeric type?
                    try:
                        val = LabelParser._convert_to_numeric(val)
                    except ValueError:
                        pass
                    # Should this value have units?
                    if match.group('units'):
                        val = Quantity(val, match.group('units'))
                    # Add it to the dict
                    obj.update({key: val})
            except StopIteration:
                break
        return obj

    @staticmethod
    def _convert_to_numeric(s):
        """Convert a string to its appropriate numeric type"""
        if re.match(INTEGER, s):
            return int(s)
        elif re.match(FLOATING, s):
            return float(s)
        else:
            raise ValueError


#-------


def menu_import(self, context):
    self.layout.operator(ImportHiRISETerrain.bl_idname, text=ImportHiRISETerrain.bl_label)


classes = (
    ImportHiRISETerrain,
    TerrainPanel,
    ReloadTerrain,
)


def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    bpy.types.TOPBAR_MT_file_import.append(menu_import)


def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    bpy.types.TOPBAR_MT_file_import.remove(menu_import)


if __name__ == '__main__':
    register()
