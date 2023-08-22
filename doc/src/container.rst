Run precon in a container
=========================

It is possible, and sometimes beneficial, to run precon in a container. If the container is used on a single machine the license file can be mapped into the container and precon can be activated as describend above.
If the container should be distributed on multiple machines a floating license needs to be added to the container as environment variable :ref:`see above <Activate-precon-to-run-on-multiple-machines>`.

Build and run the container on a single machine
-----------------------------------------------

1. Create a Dockerfile, which installs precon and copies the recon script into it. For example:
    .. code-block:: Dockerfile

        FROM python:3.9

        RUN pip install gt-precon --extra-index-url https://precon:<ACCESS_TOKEN>@pypi.gyrotools.com/simple/
        COPY ./examples/simple_recon.py /

   
2. Build the container:
   
   .. code-block:: bash
    
        docker build -t precon .

   
3. Create a license file on the host in a location of your  choice:

    .. code-block:: bash
        
        mkdir /home/gyrotools/license
        touch /home/gyrotools/license/license.key

   
4. Map the license file into the container to `/etc/gyrotools` and activate precon:
   
    .. code-block:: bash
    
        docker run --rm -v /home/gyrotools/license:/etc/gyrotools precon python -m precon license --activate <ACTIVATION_TOKEN>


5. Run the recon by mapping the license file and the data into the container:
    
    .. code-block:: bash

        docker run --rm -v /home/gyrotools/license:/etc/gyrotools -v /home/gyrotools/data:/data precon python /simple_recon.py /data/rawfile.raw

   
Use the container on multiple machines
--------------------------------------

1. Get a floating license key as described :ref:`above <Activate-precon-to-run-on-multiple-machines>`

2. Add the key as environment variable to the container (in the build or run stage):
   
    .. code-block:: bash
        
        docker run --rm -e PRECON_LICENSE_KEY=<YOUR_LICENSE_KEY> -v /home/gyrotools/data:/data precon python /simple_recon.py /data/rawfile.raw
