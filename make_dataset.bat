:: ===========================
:: ===========================
:: ===========================

@set exe=%USERPROFILE%\Desktop\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\net5.0\dimorphics_dataset.exe
@set output=e:\dataset2
@md %output%\

:: ===========================
:: ===========================
:: ===========================

set pos_class_id=+1
set pos_class_name=dimorphic_coil

set neg_class_id=-1
set neg_class_name=standard_coil

:: ===========================
:: ===========================
:: ===========================

set min_sequence_length=3
set use_dssp3=true
set max_features=100
set use_children=true
set verbose=true

set class_id=%pos_class_id%
set class_name=%pos_class_name%

:: ===========================
:: ===========================
:: ===========================



:: =========================
:: == 1d features dataset ==
:: =========================
set dimension=1

:: ===================================
:: == 1d interface subsequence area ==
:: ===================================
set area=i
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt

:: ===========================
:: == 1d neighbourhood area ==
:: ===========================
set area=n
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt

:: ===========================
:: == 1d protein chain area ==
:: ===========================
set area=p
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt



:: =========================
:: == 2d features dataset ==
:: =========================
set dimension=2

:: ===================================
:: == 2d interface subsequence area ==
:: ===================================
set area=i
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt

:: ===========================
:: == 2d neighbourhood area ==
:: ===========================
set area=n
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt

:: ===========================
:: == 2d protein chain area ==
:: ===========================
set area=p
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt



:: =========================
:: == 3d features dataset ==
:: =========================
set dimension=3

:: ===================================
:: == 3d interface subsequence area ==
:: ===================================
set area=i
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt

:: ===========================
:: == 3d neighbourhood area ==
:: ===========================
set area=n
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt

:: ===========================
:: == 3d protein chain area ==
:: ===========================
set area=p
%exe% -area=%dimension%%area% -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> %output%\stdout_%dimension%%area%_(%class_id%)_(%class_name%).txt 2> %output%\stderr_%dimension%%area%_(%class_id%)_(%class_name%).txt




:: ===========================
:: ===========================
:: ===========================
pause
