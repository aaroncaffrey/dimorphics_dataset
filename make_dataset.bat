md e:\dataset\2i\dimorphic_coil\
md e:\dataset\2n\dimorphic_coil\
md e:\dataset\2p\dimorphic_coil\
md e:\dataset\3i\dimorphic_coil\
md e:\dataset\3n\dimorphic_coil\
md e:\dataset\3p\dimorphic_coil\
md e:\dataset\2i\standard_coil\
md e:\dataset\2n\standard_coil\
md e:\dataset\2p\standard_coil\
md e:\dataset\3i\standard_coil\
md e:\dataset\3n\standard_coil\
md e:\dataset\3p\standard_coil\


:: 2d

:: DC 2d interface subsequence area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=2i -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\2i\dimorphic_coil\ 1> e:\dataset\2i\dimorphic_coil\stdout.txt 2> e:\dataset\2i\dimorphic_coil\stderr.txt
:: SC 2d interface subsequence area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=2i -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\2i\standard_coil\ 1> e:\dataset\2i\standard_coil\stdout.txt 2> e:\dataset\2i\standard_coil\stderr.txt


:: DC 2d neighbourhood area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=2n -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\2n\dimorphic_coil\ 1> e:\dataset\2n\dimorphic_coil\stdout.txt 2> e:\dataset\2n\dimorphic_coil\stderr.txt
:: SC 2d neighbourhood area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=2n -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\2n\standard_coil\ 1> e:\dataset\2n\standard_coil\stdout.txt 2> e:\dataset\2n\standard_coil\stderr.txt


:: DC 2d protein area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=2p -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\2p\dimorphic_coil\ 1> e:\dataset\2p\dimorphic_coil\stdout.txt 2> e:\dataset\2p\dimorphic_coil\stderr.txt
:: SC 2d protein area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=2p -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\2p\standard_coil\ 1> e:\dataset\2p\standard_coil\stdout.txt 2> e:\dataset\2p\standard_coil\stderr.txt


:: 3d

:: DC 3d interface subsequence area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=3i -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\3i\dimorphic_coil\ 1> e:\dataset\3i\dimorphic_coil\stdout.txt 2> e:\dataset\3i\dimorphic_coil\stderr.txt
:: SC 3d interface subsequence area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=3i -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\3i\standard_coil\ 1> e:\dataset\3i\standard_coil\stdout.txt 2> e:\dataset\3i\standard_coil\stderr.txt


:: DC 3d neighbourhood area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=3n -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\3n\dimorphic_coil\ 1> e:\dataset\3n\dimorphic_coil\stdout.txt 2> e:\dataset\3n\dimorphic_coil\stderr.txt
:: SC 3d neighbourhood area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=3n -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\3n\standard_coil\ 1> e:\dataset\3n\standard_coil\stdout.txt 2> e:\dataset\3n\standard_coil\stderr.txt


:: DC 3d protein area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=3p -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\3p\dimorphic_coil\ 1> e:\dataset\3p\dimorphic_coil\stdout.txt 2> e:\dataset\3p\dimorphic_coil\stderr.txt
:: SC 3d protein area:
"C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe" -area=3p -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=e:\dataset\3p\standard_coil\ 1> e:\dataset\3p\standard_coil\stdout.txt 2> e:\dataset\3p\standard_coil\stderr.txt


pause
