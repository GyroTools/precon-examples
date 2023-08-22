Getting Started
===============

Import precon:

.. code-block:: python

    import precon as pr


Read the parameters from a raw file (the data must have been acquired with the ReconFrame patch):

.. code-block:: python
    
    pars = pr.Parameter(Path('my_rawfile.raw'))


Obtain a parameter:

.. code-block:: python
    
    vs = pars.get_voxel_sizes(mix=0)


Obtain a GOAL parameter by name:

.. code-block:: python
    
    flip_angle = pars.get_value('EX_ACQ_flip_angle', default=[90.0])


Please note that `get_value` always returns a list of values. If the parameter is not found then the specified default value is returned. If no default value is given `None` is returned.

Obtain a GOAL object attribute:

.. code-block:: python

    sq_base = pars.goal.get_object('SQ', 'base')
    sq_base_dur = sq_base.get_attribute('dur', cmp=0)


Read the raw file and create k-space data:

.. code-block:: python

    parameter2read = pr.Parameter2Read(pars.labels)
    for mix in parameter2read.mix:
        for stack in parameter2read.stack:
            parameter2read.stack = stack
            parameter2read.mix = mix

            # read data
            with open(args.rawfile, 'rb') as raw:
                data, labels = pr.read(raw, parameter2read, pars.labels, pars.coil_info)

            # sort and zero fill data (create k-space)
            cur_recon_resolution = pars.get_recon_resolution(mix=mix, xovs=False, yovs=True, zovs=True)
            data, labels = pr.sort(data, labels, output_size=cur_recon_resolution)


Please note that the read function only reads data of the same size and with the same geometry. Therefore, you should always loop over the number of mixes and stacks. 

Examples of a complete reconstruction can be found `here <https://github.com/GyroTools/precon-examples/tree/master/examples>`_.

Example data can be found `here <https://github.com/GyroTools/precon-examples/tree/master/data>`_.


