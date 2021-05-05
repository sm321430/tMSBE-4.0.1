
import shutil
import sys
import glob

input1 = sys.argv

number_k_points = 100
if len(input1)>1:
	number_k_points = input1[1]


if len(glob.glob('dynamicInterpolationTable_e_in_MASTER_'+str(number_k_points)))==0 or len(glob.glob('dynamicInterpolationTable_e_out_MASTER_'+str(number_k_points)))==0 or len(glob.glob('dynamicInterpolationTable_h_in_MASTER_'+str(number_k_points)))==0 or len(glob.glob('dynamicInterpolationTable_h_out_MASTER_'+str(number_k_points)))==0:
	print 'Cannot find the files: dynamicInterpolationTable_e_in_MASTER_'+str(number_k_points)
	print 'Cannot find the files: dynamicInterpolationTable_e_out_MASTER_'+str(number_k_points)
	print 'Cannot find the files: dynamicInterpolationTable_h_in_MASTER_'+str(number_k_points)
	print 'Cannot find the files: dynamicInterpolationTable_h_out_MASTER_'+str(number_k_points)
	sys.exit()
else:
	print 'Duplicating MASTER_' + str(number_k_points) +' files '



for i in range(1,11):
	baseString = 'dynamicInterpolationTable_e_in_QW'
	newFile = baseString + str(i) + '.dat'
	shutil.copy2('dynamicInterpolationTable_e_in_MASTER_'+str(number_k_points), newFile)
	
	baseString = 'dynamicInterpolationTable_e_out_QW'
	newFile = baseString + str(i) + '.dat'
	shutil.copy2('dynamicInterpolationTable_e_out_MASTER_'+str(number_k_points), newFile)
	
	baseString = 'dynamicInterpolationTable_h_in_QW'
	newFile = baseString + str(i) + '.dat'
	shutil.copy2('dynamicInterpolationTable_h_in_MASTER_'+str(number_k_points), newFile)
	
	baseString = 'dynamicInterpolationTable_h_out_QW'
	newFile = baseString + str(i) + '.dat'
	shutil.copy2('dynamicInterpolationTable_h_out_MASTER_'+str(number_k_points), newFile)
