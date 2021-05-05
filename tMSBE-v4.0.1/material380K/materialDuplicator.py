
import shutil



baseString = 'material_QW'

for i in range(2,30):
	newFile = baseString + str(i) + '.config'
	shutil.copy2('material_QW1.config', newFile)
