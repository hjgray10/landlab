import numpy as np

from landlab import CLOSED_BOUNDARY, INACTIVE_LINK, Component, FieldError
from landlab.grid.divergence import calc_flux_div_at_node


class BagnoldWindTransporter(Component):
    """
    Component to simulate the movement of sediment by wind action.

    The component uses the Bagnold 1937 formula to calculate the changes in
    the field stored at "soil__depth". Here soil refers generally to material
    that can be entrained by wind.

    TODO: I would put a statement here about expectations what fields this
    component uses as well as a statement about wind direction units and
    reference frame.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid.
    parameter_name : type
        Description (including default value)

        TODO: add all parameters.

    Examples
    --------
    >>> from landlab.components import BagnoldWindTransporter
    >>> # eventually we will put an example here.
    """

    _name = "BagnoldWindTransporter"

    _input_var_names = ("soil__depth", "wind__direction" "wind__shear_velocity")

    _output_var_names = ("soil__depth", "wind_sed__flux")

    _var_units = {
        "soil__depth": "m",
        "wind__direction": "TODO"
        "wind_sed__flux": "m3/s",
        "wind__shear_velocity": "m/s",
    }

    _var_mapping = {
        "soil__depth": "node",
        "wind__direction": "node"
        "wind__shear_velocity": "link",
        "wind__sed_flux": "link",
    }

    _var_doc = {
        "soil__depth": "depth of soil/weather bedrock",
        "wind__sed_flux": "local wind-blown sediment flux from Bagnold 1937",
        "wind__shear_velocity": "the shear velocity at a node",
        "wind__direction": "TODO: state reference frame here too"
    }

        
    def __init__(self, 
                 grid,
        grain_size=0.001,  # put in default values
        sediment_density=3000,  # default values.
        **kwds
    ):
        super(BagnoldWindTransporter, self).__init__(grid)
        
        self._grid = grid

        # Unless it is really important that a user acesses the grain size
        # I would recommend saving it as a private variable (put the underscore first)
        #  also, I don't think you use the grain size with this as written.
        self._d = grain_size
        self._rho = sediment_density
    
        # Check for soil depth field and build array if necessary
        if "soil__depth" in grid.at_node:
            self.soil__depth = grid.at_node["soil__depth"]
        else:
            self.soil__depth = grid.add_zeros("soil__depth", at="node", dtype=float)
            
        # Create a wind sediment flux field and build array if necessary
        if "wind__sed_flux" in self.grid.at_link:
            self.flux = self.grid.at_link["wind__sed_flux"]
        else:
            self.flux = self.grid.add_zeros("link", "wind__sed_flux")
            
        # Check for a shear velocity field and throw an error if one isn't present
        if 'wind__shear_velocity' in grid.at_node:
            self.u_star = grid.at_node['shear__velocity']
        else:
            raise FieldError(
                    'You must supply a grid field of wind shear velocities at nodes to run this component')
        
        # Check for a wind direction field
        if 'wind__direction' in grid.at_node:
            self.wind_dir = grid.at_node['wind__direction']
        else:
            raise FieldError(
                    'You must supply a grid field of wind directions at nodes to run this component')
        
    
    
    def convert_azimuth_to_cartesian(self, grid, wind_direction):
        #### DANGER DANGER. I would mot re-write the wind direction field. I would either:
        ## a) expect the wind direction in the format you want it. The end.
        # or
        # b) ask for the wind direction in a specific format and then convert it
        # to a more usable form but store that NOT in a field but in a private variable (self._wind_dir_math?)

        # note also that I've created a PR to create a function that will do this conversion for everyone.
        ## probably want to convert into cartesian math degree form
        self.grid['node']['wind__direction'] = np.abs(np.abs(self.grid['node']['wind__direction'] -360) + 90)
        self.grid['node']['wind__direction'][np.where(self.grid['node']['wind__direction']>359)] -= 360
        
        ## then convert into radians
        self.grid['node']['wind__direction'] = np.deg2rad(self.grid['node']['wind__direction'])
    


    def wind_sed_flux_calc(self, grid, grain_size, sediment_density):
        
        ## Let's first break the wind vector (direction and magnitude) on the nodes into x and y components and sum them on the links
        
        # since you wrote self.u_star = self.grid.at_node['wind__shear_velocity'] earlier, you can now use that instead of
        # self.grid['node']['wind__shear_velocity']
                axis=1)/2.0
                
        # ditto for wind__direction and the at link fields at the end of this function.
                self.grid['node']['wind__shear_velocity'][self.grid.nodes_at_link]*np.sin(self.grid['node']['wind__direction'][self.grid.nodes_at_link]),
                axis=1)/2.0
               
        ## Now let's figure out what the angle of the wind is on the links         
                
        arctan_of_wind_x_y_components = np.arctan(wind_avg_y_component_on_link/wind_avg_x_component_on_link)
        
        ## figure out which math quadrant the wind is pointing in
        
        quadrant_of_wind_arctan = np.zeros(np.shape(wind_avg_x_component_on_link))
        
        # x = positive, y = positive : first quadrant
        quadrant_of_wind_arctan[np.logical_and(wind_avg_x_component_on_link>0, wind_avg_y_component_on_link>0)] = 1
        
        # x = negative, y = positive : second quadrant
        quadrant_of_wind_arctan[np.logical_and(wind_avg_x_component_on_link<0, wind_avg_y_component_on_link>0)] = 2
        
        # x = negative, y = negative : third quadrant
        quadrant_of_wind_arctan[np.logical_and(wind_avg_x_component_on_link<0, wind_avg_y_component_on_link<0)] = 3
        
        # x = positive, y = negative : fourth quadrant
        quadrant_of_wind_arctan[np.logical_and(wind_avg_x_component_on_link>0, wind_avg_y_component_on_link<0)] = 4
        
        # convert arctan result to the wind angle in math orientation
        wind_angle_on_the_links = arctan_of_wind_x_y_components
        wind_angle_on_the_links[quadrant_of_wind_arctan==1] += 0.0
        wind_angle_on_the_links[quadrant_of_wind_arctan==2] += np.pi #180.0
        wind_angle_on_the_links[quadrant_of_wind_arctan==3] += np.pi #180.0
        wind_angle_on_the_links[quadrant_of_wind_arctan==4] += 2.0*np.pi #360.0
        
        ## Great, now lets get the magitude of the wind on the links and calculate the sed flux
        
        wind_shear_on_the_links = np.sqrt(wind_avg_x_component_on_link**2.0 + wind_avg_y_component_on_link**2.0)
        
        wind_sed_flux_mag_on_link = (1.0 / self._rho) * (
            1.8
            * (1.225 / (9.81))
            * np.sqrt(self._d / (2.50e-6))
            * np.power(wind_shear_on_the_links, 3)
        )
        
        angle_factor = np.cos(wind_angle_on_the_links - self.grid.angle_of_link)
    #    angle_factor[angle_factor==np.cos(np.pi/2.0)] = 0.0
        wind_sed_flux_on_link = wind_sed_flux_mag_on_link*angle_factor
        wind_sed_flux_on_link[
            self.grid.status_at_link == INACTIVE_LINK
        ] = (
            0
        )  ## DANGER. Don't use the number 4, import the correct variable from Landlab.
        self.grid['link']['wind__link_shear_velocity'] = wind_shear_on_the_links
        self.grid['link']['wind__angle_on_the_link'] = wind_angle_on_the_links
        self.grid['link']['wind__sed_mass_flux'] = wind_sed_flux_on_link 
            
        # danger danger. you want to use [:] at the end of this.
        # eg.  self.grid['link']['wind__link_shear_velocity'][:] = wind_shear_on_the_links
        # that way you replace the contents of the array without changing the pointer.
        
    
    ## Ok now need to start thinking about boundary conditions
    
    def zero_soil_nodes(self, grid, dt):
        
        ## Need the angle of the link from the perpective of the node
        
        incoming_link_IDs = self.grid.links_at_node[self.grid.link_dirs_at_node==1]
        outgoing_link_IDs =  self.grid.links_at_node[self.grid.link_dirs_at_node==-1]
        
        incoming_link_angles = self.grid.angle_of_link_about_head[incoming_link_IDs]
        outgoing_link_angles = self.grid.angle_of_link[outgoing_link_IDs]
        
        link_angle_from_nodes_perspective = np.nan*np.zeros(np.shape(self.grid.links_at_node))
        link_angle_from_nodes_perspective[self.grid.link_dirs_at_node==1] = incoming_link_angles
        link_angle_from_nodes_perspective[self.grid.link_dirs_at_node==-1] = outgoing_link_angles
        
        wind_angles_reshape = np.nan*np.zeros(np.shape(self.grid.links_at_node))
        wind_angles_reshape[self.grid.links_at_node!=-1] = self.grid['link']['wind__angle_on_the_link'][self.grid.links_at_node[self.grid.links_at_node!=-1]]
        
        wind_incoming_or_outgoing = np.round(np.cos(wind_angles_reshape - link_angle_from_nodes_perspective)) 
        # 1 is outgoing, 0 is no flow, -1 is incoming, nan is an outside link
        
        # create a boolean array for incoming and outgoing links and close closed boundary links
        incoming_links_bool = np.logical_and(
            np.array([wind_incoming_or_outgoing == -1]),
            self.grid.link_status_at_node != CLOSED_BOUNDARY,
        )  # np.array([wind_incoming_or_outgoing==-1]) #
        outgoing_links_bool = np.logical_and(
            np.array([wind_incoming_or_outgoing == 1]),
            self.grid.link_status_at_node != CLOSED_BOUNDARY,
        
        # get the nodes that will go negative next timestep - "danger nodes"
        
        ##### note: need to fix for when np.max results in a zero #######
        
        outgoing_flux_per_node = self.grid['link']['wind__sed_mass_flux'][self.grid.links_at_node] # get the sed flux at each link w.r.t. the node
        outgoing_flux_per_node[incoming_links_bool[0,:,:]] = 0.0 # for this array, set the incoming flux to zero as we don't care about it
        total_outgoing_flux_per_node = np.sum(np.abs(outgoing_flux_per_node), axis=1) # figure out the outgoing flux
        
        zero_soil_danger_nodes = np.where(np.logical_and(self.grid['node']['soil__depth'] <= total_outgoing_flux_per_node*dt, self.grid.status_at_node!=4)) # these are the danger nodes
        
        # pick out the outgoing links on the danger nodes
        
        flux_on_links_at_danger_nodes = outgoing_flux_per_node[zero_soil_danger_nodes,:]
        
        # redistribute the sediment on the node into the downwind links relative to the outgoing flux
        
        relative_flux_on_danger_links_at_danger_nodes = (
                flux_on_links_at_danger_nodes / np.max(np.abs(flux_on_links_at_danger_nodes[0,:,:]),axis=1)[:,None])
        
        soil__depth_in_the_danger_nodes = self.grid['node']['soil__depth'][zero_soil_danger_nodes]
        new_outgoing_flux_on_links_out_of_danger_nodes = relative_flux_on_danger_links_at_danger_nodes*(soil__depth_in_the_danger_nodes/dt)[:,None]
        
        # now need to put it back into the larger flux array without screwing up the incoming fluxes
        
        replacement_link_index = self.grid.links_at_node[zero_soil_danger_nodes,][outgoing_links_bool[0,zero_soil_danger_nodes,]]
        replacement_link_flux = new_outgoing_flux_on_links_out_of_danger_nodes[outgoing_links_bool[0,zero_soil_danger_nodes,]]
        
        self.grid['link']['wind__sed_mass_flux'][replacement_link_index] = replacement_link_flux
        
    
    def bagnold_wind_transporter(self, grid, dt):
        
        dhdt = -self.grid.calc_flux_div_at_node(self.grid['link']['wind__sed_mass_flux'])
        
        self.grid['node']['soil__depth'] += dhdt*dt
         
        self.grid['node']['soil__depth'][self.grid['node']['soil__depth']<0] = 0.0
    
    
    
    def run_one_step(self,grid, grain_size, sediment_density, dt,**kwds):
        
        self.convert_azimuth_to_cartesian()
        self.windsedflux(dt)
        self.zero_soil_nodes()
        self.bagnold_wind_transporter()



































