
import glob
import sys
import math


def file_get_number_rows_cols(fname):
	num_col = 0
	num_row = 0
	with open(fname) as f:
		for num_row, l in enumerate(f):
			num_col = len(l.split());
			pass
		# if file is empty
		if (num_col > 0):
			num_row = num_row+1 # add one to get correct number of lines
	return num_row,num_col
	
def file_compare_number_rows_cols(fname, num_row, num_col):
	count = 0
	with open(fname) as f:
		for line in f:
			count = count + 1;
			num = len(line.split());
			if (num != num_col):
				print 'file: ' + fname + ' does not have a consistent number of columns.'
				print 'First noticed on line ' + str(count)
				print 'Expecting number of columns: ' + str(num_col)
				print 'Found number of columns: ' + str(num)
				sys.exit()
	if count != num_row:
		print 'file: ' + fname + ' does not have a consistent number of rows: '
		print 'number of rows found = ' + str(count)
		print 'number of rows expected = ' + str(num_row)
		sys.exit()
	#else:
		#print 'file: ' + fname + ' has size ' + str(count) + ' x ' + str(num_col)
		
def file_does_point_exist_in_master(point,master_name):
	TOL = 1.0e-12
	new_point = 1
	with open(master_name) as g:
		for line_g in g:
			row = line_g.split();
			row_id = row[0:3]
			#diff = point - row_id
			diff = [float(i) - float(j) for i, j in zip(point, row_id)];
			#print point
			#print row_id
			#print diff
			dist = math.sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2])
			if (dist < TOL):
				new_point = 0
				pass
	return new_point

print 'interpolationCheck: Start check of dynamicInterpolationTable files'

# Load all files from directory
dir_e_in  = glob.glob('dynamicInterpolationTable_e_in_*.dat');
dir_e_out = glob.glob('dynamicInterpolationTable_e_out_*.dat');
dir_h_in  = glob.glob('dynamicInterpolationTable_h_in_*.dat');
dir_h_out = glob.glob('dynamicInterpolationTable_h_out_*.dat');
num_e_in   = len(dir_e_in);
num_e_out  = len(dir_e_out);
num_h_in   = len(dir_h_in);
num_h_out  = len(dir_h_out);

# Extract all of device names
device_names = [name.replace("dynamicInterpolationTable_e_in_","").replace(".dat","") for name in dir_e_in]

#===================================================================================
# Check that the number of files loaded matches the number of files
if (num_e_in != num_e_out) or (num_h_in != num_h_out) or (num_e_in != num_h_in):
	print 'interpolationCheck: Missing or extra device files!'
	print '# e_in  = ' + str(num_e_in)
	print '# e_out = ' + str(num_e_out)
	print '# h_in  = ' + str(num_h_in)
	print '# h_out = ' + str(num_h_out)
	print 'Quitting!'
	sys.exit()
else:
	print 'interpolationCheck: Total number of files ok'
	



#===================================================================================
# Check that the grouping of files per device is consistent
for device in device_names:
	
	all_files = glob.glob('dynamicInterpolationTable_*_' + device + '.dat')
	if len(all_files) != 4:
		print 'interpolationCheck: Missing datafiles!'
		print 'Cannot find all datafiles for device: ' + device
		for el in all_files:
			print 'found file: ' + el
		sys.exit()
		
	
print 'interpolationCheck: Number of files per device ok'
#===================================================================================
# Check that the number of points in each device matches
device_number_rows = []
device_number_cols = []
for device in device_names:
	# Open one file to check number of rows and columns, this has to match for ALL files in set
	fname = 'dynamicInterpolationTable_e_in_' + device + '.dat'
	num_col = 0
	num_row = 0
	num_row, num_col = file_get_number_rows_cols(fname)
	
	# Check that each file has the correct number of rows and columns
	fname = 'dynamicInterpolationTable_e_in_' + device + '.dat'
	file_compare_number_rows_cols(fname, num_row, num_col)
		
	fname = 'dynamicInterpolationTable_e_out_' + device + '.dat'
	file_compare_number_rows_cols(fname, num_row, num_col)
		
	fname = 'dynamicInterpolationTable_h_in_' + device + '.dat'
	file_compare_number_rows_cols(fname, num_row, num_col)
		
	fname = 'dynamicInterpolationTable_h_out_' + device + '.dat'
	file_compare_number_rows_cols(fname, num_row, num_col)
	
	device_number_rows.append(num_row)
	device_number_cols.append(num_col)
	print 'interpolationCheck: Device (' + device + ') files have consistent number of rows and columns: '

#===================================================================================
# Check that each point that appears in a file is also duplicated to the other files
for device in device_names:
	# Open one file to check number of rows and columns, this has to match for ALL files in set
	fname_e_in = 'dynamicInterpolationTable_e_in_' + device + '.dat'
	fname_e_out = 'dynamicInterpolationTable_e_out_' + device + '.dat'
	fname_h_in = 'dynamicInterpolationTable_h_in_' + device + '.dat'
	fname_h_out = 'dynamicInterpolationTable_h_out_' + device + '.dat'
	
	# All file have the same size, no problem opening at the same time
	with open(fname_e_in) as f_e_in,open(fname_e_out) as f_e_out,open(fname_h_in) as f_h_in,open(fname_h_out) as f_h_out:
		for line_e_in,line_e_out,line_h_in,line_h_out in zip(f_e_in,f_e_out,f_h_in,f_h_out):
			
			row_e_in  = line_e_in.split()[0:3]
			row_e_out = line_e_out.split()[0:3]
			row_h_in  = line_h_in.split()[0:3]
			row_h_out = line_h_out.split()[0:3]
			
			diff_e  = sum([abs(float(i) - float(j)) for i, j in zip(row_e_in, row_e_out)]);
			diff_h  = sum([abs(float(i) - float(j)) for i, j in zip(row_h_in, row_h_out)]);
			diff_eh = sum([abs(float(i) - float(j)) for i, j in zip(row_e_in, row_h_in)]);
			
			if diff_e != 0 or diff_h != 0 or diff_eh != 0:
				print 'points in files dont match'
				print '(' + device + '): e_in  = ' + row_e_in
				print '(' + device + '): e_out = ' + row_e_out
				print '(' + device + '): h_in  = ' + row_h_in
				print '(' + device + '): h_out = ' + row_h_out
				print '(' + device + '): diff_e  = ' + str(diff_e)
				print '(' + device + '): diff_h  = ' + str(diff_h)
				print '(' + device + '): diff_eh = ' + str(diff_eh)
				sys.exit()
						
	print 'interpolationCheck: Device (' + device + ') has the same points in all files.'
	
#===================================================================================
#print device_names
#print device_number_rows
#print device_number_cols
device_number_cols_unique = list(set(device_number_cols)) # Recover unique number of columns
#===================================================================================
# Check if MASTER files exist and create for each number of momentum values
for data in device_number_cols_unique:
	
	fname = 'dynamicInterpolationTable_e_in_MASTER_' + str(data-3)
	dir_e_in_M  = glob.glob(fname);
	if len(dir_e_in_M)==0:
		print 'File: ' + fname + ' not found, creating blank file'
		f = open(fname,'w+')
		f.close()
		
	fname = 'dynamicInterpolationTable_e_out_MASTER_' + str(data-3)
	dir_e_out_M  = glob.glob(fname);
	if len(dir_e_out_M)==0:
		print 'File: ' + fname + ' not found, creating blank file'
		f = open(fname,'w+')
		f.close()
		
	fname = 'dynamicInterpolationTable_h_in_MASTER_' + str(data-3)
	dir_h_in_M  = glob.glob(fname);
	if len(dir_h_in_M)==0:
		print 'File: ' + fname + ' not found, creating blank file'
		f = open(fname,'w+')
		f.close()
		
	fname = 'dynamicInterpolationTable_h_out_MASTER_' + str(data-3)
	dir_h_out_M = glob.glob(fname);
	if len(dir_h_out_M)==0:
		print 'File: ' + fname + ' not found, creating blank file'
		f = open(fname,'w+')
		f.close()
	
	
	
	
	# Open one file to check number of rows and columns, this has to match for ALL files in set
	fname = 'dynamicInterpolationTable_e_in_MASTER_' + str(data-3)
	num_col = 0
	num_row = 0
	num_row, num_col = file_get_number_rows_cols(fname)
	
	
	# Check that each file has the correct number of rows and columns
	fname = 'dynamicInterpolationTable_e_in_MASTER_' + str(data-3)
	file_compare_number_rows_cols(fname, num_row, num_col)
		
	fname = 'dynamicInterpolationTable_e_out_MASTER_' + str(data-3)
	file_compare_number_rows_cols(fname, num_row, num_col)
		
	fname = 'dynamicInterpolationTable_h_in_MASTER_' + str(data-3)
	file_compare_number_rows_cols(fname, num_row, num_col)
		
	fname = 'dynamicInterpolationTable_h_out_MASTER_' + str(data-3)
	file_compare_number_rows_cols(fname, num_row, num_col)
	
	print 'interpolationCheck: MASTER_' + str(data-3) + ' files have consistent number of rows and columns'


#===================================================================================
# Merge device files with master files
for i, val in enumerate(device_names):
	
	fname_e_in 			= 'dynamicInterpolationTable_e_in_' + device_names[i] + '.dat'
	master_name_e_in 	= 'dynamicInterpolationTable_e_in_MASTER_' + str(device_number_cols[i]-3)
	fname_e_out 		= 'dynamicInterpolationTable_e_out_' + device_names[i] + '.dat'
	master_name_e_out 	= 'dynamicInterpolationTable_e_out_MASTER_' + str(device_number_cols[i]-3)
	fname_h_in 			= 'dynamicInterpolationTable_h_in_' + device_names[i] + '.dat'
	master_name_h_in 	= 'dynamicInterpolationTable_h_in_MASTER_' + str(device_number_cols[i]-3)
	fname_h_out 		= 'dynamicInterpolationTable_h_out_' + device_names[i] + '.dat'
	master_name_h_out 	= 'dynamicInterpolationTable_h_out_MASTER_' + str(device_number_cols[i]-3)
	num_new = 0
	with open(fname_e_in) as f_e_in,open(fname_e_out) as f_e_out,open(fname_h_in) as f_h_in,open(fname_h_out) as f_h_out:
		for line_e_in,line_e_out,line_h_in,line_h_out in zip(f_e_in,f_e_out,f_h_in,f_h_out):
			row = line_e_in.split();
			row_id = row[0:3]
			
			# Check if this matches anything in MASTER
			new_point = file_does_point_exist_in_master(row_id,master_name_e_in)
			
			if (new_point == 1):
				# add new point
				num_new = num_new + 1;
				with open(master_name_e_in, "a") as master_file_e_in,open(master_name_e_out, "a") as master_file_e_out,open(master_name_h_in, "a") as master_file_h_in,open(master_name_h_out, "a") as master_file_h_out:
					master_file_e_in.write(line_e_in)
					master_file_e_out.write(line_e_out)
					master_file_h_in.write(line_h_in)
					master_file_h_out.write(line_h_out)
	
	
	
	print 'interpolationCheck: Merging ' + device_names[i] + ' files to MASTER_' + str(device_number_cols[i]-3) + ' files. New points = ' + str(num_new)
					
			

























