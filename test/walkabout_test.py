import sys
import os
import subprocess
import argparse

def test_num(test_results, ver_answers):
	#testfiles format is['test results', 'answers']
	#feel free to add any doc with values in this format to be tested
	#compatible with any text based document if values are spaced from any other # character
	testfiles = [test_results, ver_answers]
	data = []
	k = 0
	#open results calculated
	while k < (len(testfiles)):
		with open(testfiles[k]) as f:
			k+=1
			results = []
			for line in f:
				line = line.split()             
		 # clean blank lines
				if line:
					for i in range(len(line)):
		 # to remove all non float values from list         
					   try:
							line = [float(i) for i in line]
					   except ValueError:
							line.remove(i)
				for x in line:
					results.append(x)
	#open answers to compare to results
		with open(testfiles[k]) as d:
			k+=1
			answers = []
			for line in d:
				line = line.split()             
		 # clean blank lines
				if line:
					for i in range(len(line)):
		 # to remove all non float values from list          
					   try:
							line = [float(i) for i in line]
					   except ValueError:
							line.remove(i)
				for x in line:
					answers.append(x)
	#start comparing results to answers with
	#tolerance percent = tp
	#tolerance for answers of '0' = +,- tz
		tp = 0.10
		tz = 0.000001
		fail = 0
		res_sum = 0
		ans_sum = 0
		
		#comparing average of middle 70% of testfiles
		#print "len ", len(results), len(answers)
		#print "starting... ", int(0.15 * len(results)), int(len(results) - 0.15 * len(results))
		#print "starting... ", int(0.15 * len(answers)), int(len(answers) - 0.15 * len(answers))
		
		for i in results[int(0.15 * len(results)) : (int(len(results) - 0.15 * len(results)))]:
				res_sum += i
		for i in answers[int(0.15 * len(answers)) : (int(len(answers) - 0.15 * len(answers)))]:
				ans_sum += i
		res_avg = res_sum / int( 0.7 * len(results))
		ans_avg = ans_sum / int( 0.7 * len(answers))
		
		#print "avg ", res_avg, ans_avg
		
		if len(results) == len(answers):
			j = 0
			while j < len(results):
				if answers[j] < 0 and tp >= 0:
					tp = -(tp)
				if answers[j] > 0:
					tp = abs(tp)
				if answers[j] == 0:
					if (results[j] > tz or results[j] < -(tz)):
						fail =1
				else:
					if results[j] > answers[j] + answers[j] * tp or results[j] < answers[j] - answers[j] * tp:
						fail = 1
				j+=1
		else:
			fail = 1
		#testing averages of middle 70% of testfiles
		#debug for upper and lower bounds
		#print "lower ", abs(ans_avg) - abs(ans_avg * tp), abs(ans_avg) + abs(ans_avg * tp)
		if fail != 0 and (abs(res_avg) > abs(ans_avg) - abs(ans_avg * tp)) and (abs(res_avg) < abs(ans_avg) + abs(ans_avg * tp)):
			fail = 0
			print "Average testing triggered with tolderance of ", tp * 100, "%...\n"
			
		if fail == 0:
			data.append('PASSED')
		else:
			data.append('FAILED')
    #prints pass/fail
    	for i in range(len(data)):
    		print "\n..... " + data[i]


#Arguments: walkabout source, threads for omp
parser = argparse.ArgumentParser(description='Tests walkabout against verified sources')
home = os.getcwd()
os.environ['OMP_NUM_THREADS'] = '1'

parser.add_argument("src", metavar="SOURCE", nargs='?', help="walkabout executable")
args = parser.parse_args()

#set defualts
defualt_walkabout = ['walkabout']
if args.src is None:
	args.src = defualt_walkabout
	

#run tests on following directories; add as many directories as needed.
test_dirs = [r'test1a', r'test1b', r'test1c', r'test1d', r'test2', r'test3', r'test4', r'test5']
for i in range(0, len(test_dirs)):
	test_files = os.path.abspath(os.path.join(test_dirs[i], 'traj.out'))
	print "\n############ Testing walkabout on " + test_dirs[i] + " ############"
	print "Test files: \n" + test_files
	open(test_files, 'a').close()
	os.remove(test_files)
	subprocess.call(args.src, cwd=test_dirs[i])
	os.chdir(test_dirs[i])
	test_num('traj.out', 'traj.out.ans')
	os.chdir(home)
	

