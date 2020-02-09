set exe=C:\Users\Aaron\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\netcoreapp3.1\dimorphics_dataset.exe
set output=e:\dataset
md %output%\

:: 2d interface subsequence area:
:: %exe% -area=2i -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_2i_(+1)_(dimorphic_coil).txt 2> %output%\stderr_2i_(+1)_(dimorphic_coil).txt
:: %exe% -area=2i -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_2i_(-1)_(standard_coil).txt 2> %output%\stderr_2i_(-1)_(standard_coil).txt

:: 2d neighbourhood area:
:: %exe% -area=2n -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_2n_(+1)_(dimorphic_coil).txt 2> %output%\stderr_2n_(+1)_(dimorphic_coil).txt
%exe% -area=2n -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_2n_(-1)_(standard_coil).txt 2> %output%\stderr_2n_(-1)_(standard_coil).txt

:: 2d protein chain area:
%exe% -area=2p -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_2p_(+1)_(dimorphic_coil).txt 2> %output%\stderr_2p_(+1)_(dimorphic_coil).txt
%exe% -area=2p -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_2p_(-1)_(standard_coil).txt 2> %output%\stderr_2p_(-1)_(standard_coil).txt

:: 3d interface subsequence area:
%exe% -area=3i -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_3i_(+1)_(dimorphic_coil).txt 2> %output%\stderr_3i_(+1)_(dimorphic_coil).txt
%exe% -area=3i -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_3i_(-1)_(standard_coil).txt 2> %output%\stderr_3i_(-1)_(standard_coil).txt

:: 3d neighbourhood area:
%exe% -area=3n -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_3n_(+1)_(dimorphic_coil).txt 2> %output%\stderr_3n_(+1)_(dimorphic_coil).txt
%exe% -area=3n -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_3n_(-1)_(standard_coil).txt 2> %output%\stderr_3n_(-1)_(standard_coil).txt

:: 3d protein chain area:
%exe% -area=3p -use_dssp3=true -class_id=+1 -class_name=dimorphic_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_3p_(+1)_(dimorphic_coil).txt 2> %output%\stderr_3p_(+1)_(dimorphic_coil).txt
%exe% -area=3p -use_dssp3=true -class_id=-1 -class_name=standard_coil -min_sequence_length=3 -max_features=100 -output_folder=%output%\ -use_children=true -verbose=true 1> %output%\stdout_3p_(-1)_(standard_coil).txt 2> %output%\stderr_3p_(-1)_(standard_coil).txt

pause
