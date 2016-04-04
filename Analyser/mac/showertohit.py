import sys,os

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem

from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root files
for x in xrange(len(sys.argv)-1):
	my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("ShowerToHit_output.root");

# Set c++ script as analysis unit
ana_unit = fmwk.ShowerToHit()

# Add analysis unit to ana_processor instance
my_proc.add_process(ana_unit)

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Run it.
my_proc.run();

print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
