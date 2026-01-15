"""
instantiate dataset
"""
from givernylocal.turbulence_dataset import *
from givernylocal.turbulence_toolkit import *
import os

script_path = os.path.dirname(__file__)
auth_token = 'cn.edu.pku.stu.zza-5ea337ea'
dataset_title = 'isotropic1024coarse'
output_path = script_path + '/.giverny_output'

# instantiate the dataset.
dataset = turb_dataset(dataset_title = dataset_title, output_path = output_path, auth_token = auth_token)

"""
process getCutout data
"""
variable = 'velocity'


x_range = [1, 1024]
y_range = [1, 1024]
z_range = [1024, 1024]
t_range = [1, 1]

x_stride = 1
y_stride = 1
z_stride = 1
t_stride = 1

"""
use the tools and processing gizmos.
"""
# combine all of the axis data together for simplicity.
axes_ranges = np.array([x_range, y_range, z_range, t_range])
strides = np.array([x_stride, y_stride, z_stride, t_stride])

# process a brick cutout.
result = getCutout(dataset, variable, axes_ranges, strides)


"""
write the cutout results to hdf5 and xmf files
"""
output_filename = dataset_title + "_z1024"

"""
use the hdf5 and xmf writing gizmo.
"""
# writes the output hdf5 and xmf files.
write_cutout_hdf5_and_xmf_files(dataset, result, output_filename)