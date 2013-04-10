from glob import glob

ccflags=["-Wall", "-Werror", "-DNDEBUG", "-O3"]
#ccflags=["-Wall", "-Werror", "-O3"]

import os
uname = os.uname()
if uname[0] == 'Darwin' and uname[4] == 'i386':
        ccflags.append('-mdynamic-no-pic')

source_files = glob('*.c') + glob('primorials/*.c')

cpppath=['..']

StaticLibrary(target='sspar',
              source=source_files,
              CCFLAGS=ccflags,
              CPPPATH=cpppath)
