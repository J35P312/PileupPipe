import os
import subprocess
import sys

programDirectory = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(programDirectory,"config_template/FindSV.config"), 'r') as myfile:
    template=myfile.read()
        

template=template.replace("{working_dir}", programDirectory)
os.system("chmod +x {}/internal_scripts/*".format(programDirectory))


print template
