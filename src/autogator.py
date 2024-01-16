# Global Imports
import numpy as np
import FlowCal
from matplotlib import pyplot as plt

# Local Imports
from .gating_tools import scale_polygon, gate_contour


# ----------------------------------------
# Main
# ----------------------------------------

class AutoGator:
    def __init__(
        self,
        gate_fraction: float,
        negative_gate_fraction: 0.98,
        fmos= False,
        second_dimension= 'SSC-A',
        scale_factor = False,
        control_sample = False
    ):
        """
        Args:\n 
        gate_fraction: Percentage of cells to try and keep when automatically drawing gate.\n
        negative_gate_fraction: Percentage of cells to omit based on FMO before drawing the 'positive' gate on remaining.\n
        fmos: If included, must be a list of dictionaries with this format [{name: {name of omitted flourescent}, location: {location of the *.fcs file}}].  The dictionary must also contain channel_name if this is not labelled in the fcs file\n
        second_dimension: The secondary axis to use when automatically gating.  Best to use FSC-A or SSC-A\n
        scale_factor: Factor by which to increase the size of automatically generated polygon gate.  Not recommended but may be useful if minimizing Type-II errors is more important than Type-I\n
        control_sample: Location of your control sample, which should be a sample containing all all fluorescent dyes and all populations of interest. If you have multiple FMOs, this is optional but it is ideal to include one\n
        
        \nAssumptions:\n
        That labelling is consistent between batches of *.fcs files (If a channel is named 'BV786-A' in one file it cannot be labelled 'BrilliantViolet 786-A' in another)\n
        That all measurements and channels are in the same order for each FMO.  This is the weakest assumption and this should be made robust against it as soon as possible, TODO potentially by sorting by channels
        """
        self.gate_fraction = gate_fraction
        self.negative_gate_fraction = negative_gate_fraction
        self.second_dimension = second_dimension
        self.scale_factor = scale_factor
        self.fmos = {} # FlowCal sample objects for each fmo
        self.channel_lookup = {} # Usable channel name for each channel label
        if fmos:
            for fmo in fmos:
                if 'channel_name' not in fmo:
                    fmo['channel_name'] = False
                self.add_fmo(
                    name= fmo['name'],
                    location= fmo['location'],
                    channel_name= fmo['channel_name'],
                    recalculate_controls=False
                    )
            self._calculate_controls()
        if control_sample:
            self.control = FlowCal.io.FCSData(control_sample)
            self.user_supplied_control = True
        else:
            self.control = self._synthesize_control()
            self.user_supplied_control = False
    
    def _extract_channel_dict(
        self,
        sample: FlowCal.io.FCSData
    ):
        channel_labels = sample.channel_labels()
        channel_dict = {channel_labels[x]: sample.channels[x] for x in range(len(sample.channels)) if channel_labels[x]}
        for label, name in channel_dict.items():
            self.channel_lookup[label] = name
    
    def _synthesize_control(
       self
    ): 
        control = list(self.fmos).values()[0] # Hopefully the first fmo is representative
        # __array_wrap__ is a FlowCal method which presumably wraps a numpy array returning the FCSData object
        #       This seems to apply all appropriate attributes from the source object, so its a good way to concatenate the data (I hope?)
        
        control = control.__array_wrap__(np.vstack(tuple(self.fmos.values())))
        self.control = control
    
    def recalculate_gates(
        self
    ):
        for name, fmo in self.fmos.items():
            channels = [self.second_dimension, self.channel_lookup[name]]
            negatively_gated= FlowCal.gate.density2d(
                fmo,
                channels = channels,
                gate_fraction= self.negative_gate_fraction,
                full_output=True
            )
            negative_contour= negatively_gated.contour
            
            negative_mask= gate_contour(
                flowcal_object= self.control,
                channels= channels,
                contour= negative_contour
            )
            
            trimmed_control= self.control[~negative_mask]
            
            positively_gated= FlowCal.gate.density2d(
                trimmed_control,
                channels = channels,
                gate_fraction= self.gate_fraction,
                full_output=True
            )
            
            positive_contour = positively_gated.contour
            
            if self.scale_factor:
                positive_contour[0] = scale_polygon(
                    polygon= positive_contour[0],
                    scale_factor= self.scale_factor
                )
            
            # Its easier to just add attributes to the FCSData object than hash it separately
            fmo.negative_contour= negative_contour
            fmo.positive_contour= positive_contour
            
            
    
    def add_fmo(
        self,
        name: str,
        location: str,
        channel_name = False,
        recalculate_controls = True
    ):
        """
        Args:
        name: Name of omitted fluorescent, ie if dye corresponding to 'CD4' was omitted call this 'CD4'\n
        location: Location of the *.fcs file of the FMO\n
        channel_name: Name of the channel corresponding to the omitted fluorescent.  This is optional is the *.fcs file contains appropriately labelled channels\n
        recalculate_control: This can be disabled if calling this method multiple times to prevent redundant calculation when adding batches of FMOs.  Make sure to enable it for the final FMO\n
        """
        
        if channel_name:
            self.channel_lookup[name] = channel_name
    
        sample = FlowCal.io.FCSData(location)
        self._extract_channel_dict(sample)
        self.fmos[name] = sample
        
        if recalculate_controls:
            if len(self.fmos)> 1:
                self._synthesize_control()
            self.recalculate_gates()
            
        return self
    
    def gate(
        self,
        sample: FlowCal.io.FCSData,
        target: str,
        plot = False
    ) -> np.array:
        """
        Args:
        sample: FCSData() object containing the sample to be gated\n
        target: Channel label to gate for
        """
        
        positive_contour = self.fmos[target].positive_contour
        main_channel = self.channel_lookup[target]
        channels = [self.second_dimension, main_channel]
        
        gate_mask = gate_contour(
            flowcal_object= sample,
            channels= channels,
            contour= positive_contour
        )
        
        
        if plot:
            FlowCal.plot.density_and_hist(
                sample,
                gated_data=sample[gate_mask],
                gate_contour=positive_contour,
                density_channels= channels,
                hist_channels=[main_channel],
                density_params={
                    'mode': 'scatter'
                }
            )
            plt.tight_layout()
            plt.show()
        
        return gate_mask
    
    def get_gate_boundry(
        self,
        target: str,
    ) -> dict:
        """
        Args:
        sample: FCSData() object containing the sample to be gated\n
        target: Channel label whose contour to grab
        """
        positive_contour = self.fmos[target].positive_contour
        main_channel = self.channel_lookup[target]
        channels = [self.second_dimension, main_channel]
        
        return {
            'contour': positive_contour,
            'channels': channels
        }