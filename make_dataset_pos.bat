:: ===========================
:: ===========================
:: ===========================

setlocal EnableDelayedExpansion
@set exe=c:\phd\phd_work\dimorphics_dataset\dimorphics_dataset\bin\x64\Release\net6.0\dimorphics_dataset.exe
@set output=C:\phd\_march_2022_dataset\aaindex_only\pos\
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
set max_features=112
set use_children=true
set verbose=true

::set class_id=%neg_class_id%
::set class_name=%neg_class_name%

set class_id=%pos_class_id%
set class_name=%pos_class_name%

:: ===========================
:: ===========================
:: ===========================

:: for all features: blank or =1 enable, =0 disable

set features1i=-1i
set features1n=-1n
set features1p=-1p

set features2i=-2i
set features2n=-2n
set features2p=-2p

set features3i=-3i
set features3n=-3n
set features3p=-3p


:: for individually selected features: blank or =1 enable, =0 disable

set features1i=-1i.pse_aac=1 -1i.sable=1 -1i.blast_pssm=1 -1i.aaindex=1 -1i.geometry=1 -1i.iupred2a=1 -1i.stackdppred=1 -1i.r_peptides=1 -1i.r_protr=1
set features1n=-1n.pse_aac=1 -1n.sable=1 -1n.blast_pssm=1 -1n.aaindex=1 -1n.geometry=1 -1n.iupred2a=1 -1n.stackdppred=1 -1n.r_peptides=1 -1n.r_protr=1
set features1p=-1p.pse_aac=1 -1p.sable=1 -1p.blast_pssm=1 -1p.aaindex=1 -1p.geometry=1 -1p.iupred2a=1 -1p.stackdppred=1 -1p.r_peptides=1 -1p.r_protr=1

set features2i=-2i.mpsa=1
set features2n=-2n.mpsa=1
set features2p=-2p.mpsa=1

set features3i=-3i.pse_ssc_dssp=1 -3i.foldx=1 -3i.ring=1 -3i.sasa=1 -3i.tortuosity=1 -3i.intramolecular=1
set features3n=-3n.pse_ssc_dssp=1 -3n.foldx=1 -3n.ring=1 -3n.sasa=1 -3n.tortuosity=1 -3n.intramolecular=1
set features3p=-3p.pse_ssc_dssp=1 -3p.foldx=1 -3p.ring=1 -3p.sasa=1 -3p.tortuosity=1 -3p.intramolecular=1

::: todo: how many features are in ALL_DEEP?

:: disable dimensions/areas/features

set features1i=-1i.aaindex=1
set features1n=
set features1p=

set features2i=
set features2n=
set features2p=

set features3i=
set features3n=
set features3p=


:: ===================================
:: == 1d interface subsequence area ==
:: ===================================
if "%features1i%" NEQ "" (
	@set features=%features1i%
	@set area=1i
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)

:: ===========================
:: == 1d neighbourhood area ==
:: ===========================
if "%features1n%" NEQ "" (
	@set features=%features1n%
	@set area=1n
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)

:: ===========================
:: == 1d protein chain area ==
:: ===========================
if "%features1p%" NEQ "" (
	@set features=%features1p%
	@set area=1p
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)


:: ===================================
:: == 2d interface subsequence area ==
:: ===================================
if "%features2i%" NEQ "" (
	@set features=%features2i%
	@set area=2i
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)

:: ===========================
:: == 2d neighbourhood area ==
:: ===========================
if "%features2n%" NEQ "" (
	@set features=%features2n%
	@set area=2n
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)

:: ===========================
:: == 2d protein chain area ==
:: ===========================
if "%features2p%" NEQ "" (
	@set features=%features2p%
	@set area=2p
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)

:: ===================================
:: == 3d interface subsequence area ==
:: ===================================
if "%features3i%" NEQ "" (
	@set features=%features3i%
	@set area=3i
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)

:: ===========================
:: == 3d neighbourhood area ==
:: ===========================
if "%features3n%" NEQ "" (
	@set features=%features3n%
	@set area=3n
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)

:: ===========================
:: == 3d protein chain area ==
:: ===========================
if "%features3p%" NEQ "" (
	@set features=%features3p%
	@set area=3p
	@set stdout="%output%\stdout_!area!_(%class_id%)_(%class_name%).txt"
	@set stderr="%output%\stderr_!area!_(%class_id%)_(%class_name%).txt"
	@echo stdout: !stdout!
	@echo stderr: !stderr!
	%exe% !features! -use_dssp3=%use_dssp3% -class_id=%class_id% -class_name=%class_name% -min_sequence_length=%min_sequence_length% -max_features=%max_features% -output_folder=%output%\ -use_children=%use_children% -verbose=%verbose% 1> !stdout! 2> !stderr!
)



:: ===========================
:: ===========================
:: ===========================
pause
