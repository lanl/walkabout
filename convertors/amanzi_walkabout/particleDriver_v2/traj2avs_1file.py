# read walkabout output particle tracks from traj.out
# write all tracks into single AVS line 
# Usage: python traj2avs_1file.py
# Python 2.x
# tamiller@lanl.gov May 14 2018
#
# INPUT: traj.out (particle tracks written by walkabout)
# OUTPUT: traj_ptracks.inp AVS line file 
#
# Modified May 18 2018
# corrected line count to agree with header
# add warning if header count differs from number in file
#
###########################################

import sys
import os.path
import datetime
import numpy as np


def write_line_avs(infile,outfile):
    # AVS file specification: 
    # Header: <# of nodes> <# of lines> <# of attributes> <0> <0>
    # <node #> <node_x> <node_y> <node_z>
    # [step #] [2   line] [node_1] [node_2]
    # [00004  1  1  1  1]
    # [attribute_name, attribute_type]
    # [node #] [attribute_boolean_1] [attribute_boolean_2] ... [attribute_boolean_n]
    
    node_array = []     # Array for storing nodes
    tet_array = []      # Array for storing tets
    att_array = []      # Array for storing atts.
    
    pdata_in = open(infile, "r")
    nodes = []        
    line_id = []
    time_att = []  
    elev_att = []  
    
    attrib_count = 4
    tmax_step = 0
    num_ptracks = 0
    ptrack = 0
    hdrcount = 1
    num_steps = 0
    step_count = 0
    inodes = 0
    ilines = 0

    
    print("\nReading traj.out file...")

    for line in pdata_in:
        stream = line.strip().split()

        # read header information
        if 'walkabout' in stream:
              # print walkabout version
              print(stream)
              hdrcount = hdrcount+1
        elif ('tmax' in stream) or ('dtmax' in stream):
              tmax_step= int(stream[1]) 
              print("Max time steps " + str(tmax_step))
              hdrcount = hdrcount+1
        elif hdrcount == 3:
              num_ptracks = int(stream[0]) 
              hdrcount = -1
              print("Number of particle tracks " + str(num_ptracks))
              print("-------------------------------------------------------")

        # done with header now read lines of data
        # read each particle track and steps (x,y,z)
        else:

            if step_count == 0:
                ptrack = ptrack + 1
                num_steps = int(stream[0])
                ilines = ilines+num_steps-1
                print("Reading ptrack "+str(ptrack)+" steps: "+str(num_steps))
                step_count = step_count + 1

            else:

                # print(str(inodes)+": "+str(step_count)+" of "+str(num_steps)+"\n")

                time = stream[0]
                node_x = stream[1]
                node_y = stream[2]
                node_z = stream[3]

                line_id.append(ptrack)
                time_att.append(time)
                elev_att.append(node_z)

                nodes.append([node_x, node_y, node_z])
                # print("   Adding node {} of {}".format(step_count, num_steps))

                if int(step_count) == int(num_steps):
                    step_count = 0
                else:
                    step_count = step_count + 1

                inodes=inodes+1

    pdata_in.close()

    if num_ptracks != ptrack:
        print "Warning: Expected number of ptracks: "+str(num_ptracks)+" but read: "+str(ptrack) 

    print("Read Done.")

    print("Writing AVS file...")
    avs_out = open(outfile, "w")


    # Begin writing to AVS file
    attrib_count = 3
    avs_out.write("# traj.out to AVS {}\n".format(str(datetime.datetime.now())))
    avs_out.write("         {}         {}         {}         {}         {}\n".format(inodes, ilines, attrib_count, 0, 0))
    
    # write step points x,y,z
    for i in range(len(nodes)):
        ii = "{0:0>3}".format(i+1)
        avs_out.write("{}   {:.14E}  {:.14E}  {:.14E}\n".format(ii, float(nodes[i][0]), float(nodes[i][1]), float(nodes[i][2])))
    
    # write line connectivity
    # i starts at 0, line pairs starts at 1
    # skip first index when starting a new set of particle track lines

    ii = 0
    icount = 0
    for i in range(len(nodes)-1):
        ii = ii + 1
        if int(line_id[i]) == int(line_id[i+1]): 
            icount = icount +1
            avs_out.write( str(icount)+"  "+str(line_id[i])+" line  "+str(ii)+"  "+str(ii+1)+"\n")
       
    if icount != ilines:
        print "Error: Header "+str(ilines)+" differs from number of lines written: "+str(icount) 

    # write line attributes
    avs_out.write("00003  1  1  1\n")
    avs_out.write("imt1, integer \ntime, real \nelev, real \n")
        
    for i in range(len(nodes)):
        avs_out.write(str(i+1)+"  "+str(line_id[i])+"  " +
                      str(time_att[i]) + "  " + str(elev_att[i]) + "\n")
        

    print("Particle Tracks written to file: "+outfile)
    avs_out.close()
    pdata_in.close()



##############################################################################################################

def main(argv="None"):
    
    print("traj.out to AVS line file ")
    print("Executed on "+str(datetime.datetime.now()).split()[0])
    print("============================\n")
    

    infile = "traj.out" 
    if not os.path.isfile(infile):
        print("ERROR: "+infile+" doesn't exist!")
        print("Quitting...")
        sys.exit()
        
            
    outfile = "traj_ptracks.inp" 
    print("AVS Output : "+outfile)

    write_line_avs(infile, outfile)
            
    

if __name__ == "__main__":
    main()
