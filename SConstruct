import os

#~ #############################################################################
# insert a final space

GNU_OPT_FLAGS = '-Wall -O3 -g0 -I./libs '
GNU_DBG_FLAGS = '-Wall -O0 -g3 -I./libs '
AIX_OPT_FLAGS = '-O5 -qmaxmem=-1 -qarch=pwr5 -qcache=auto -qinline -qunroll -qipa -DAIX -I./libs -I/users/u/idmec/jhenriques/include '

#~ #############################################################################
base_files = [ 'DE_ExtCall.cpp', 'DESolver.cpp' ]

evars = { 'PATH' : os.environ['PATH'], 'HOME' : os.environ['HOME'] }

#~ #############################################################################
if 'DE-opt' in COMMAND_LINE_TARGETS:
  env = Environment( CXX = 'mpic++', CXXFLAGS = GNU_OPT_FLAGS, ENV = evars )
  env.Program('DE-opt', base_files, LIBS = ['config++'] )

elif 'DE-dbg' in COMMAND_LINE_TARGETS:
  env = Environment( CXX = 'mpic++', CXXFLAGS = GNU_DBG_FLAGS, ENV = evars )
  env.Program('DE-dbg', base_files, LIBS = ['config++'] )

#~ #############################################################################
# system AIX
#
elif 'xDE-opt' in COMMAND_LINE_TARGETS:
  env = Environment( CXX = 'mpic++', CXXFLAGS = AIX_OPT_FLAGS, ENV = evars )
  env.Program('xDE-opt', base_files, LIBS = [ File('/users/u/idmec/jhenriques/lib/libconfig++.a') ] )

elif 'xDE-1-opt' in COMMAND_LINE_TARGETS:
  env = Environment( CXX = 'xlc++', CXXFLAGS = AIX_OPT_FLAGS, ENV = evars )
  env.Program('xDE-1-opt', base_files, LIBS = [ File('/users/u/idmec/jhenriques/lib/libconfig++.a') ] )

#~ #############################################################################
else:
  print "No defined targets"
  exit(0)

print "BUILD_TARGETS is", map( str, BUILD_TARGETS )

#~ ###EOF#######################################################################


