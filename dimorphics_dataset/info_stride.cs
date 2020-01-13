using System.Collections.Generic;
using System.IO;

namespace dimorphics_dataset
{

    public static class column_value_extensions
    {
        public static string column_value(this string column_format_line, int start_position, int end_position, bool trim = true)
        {
            var ret = column_format_line.Substring(start_position - 1, (end_position - start_position) + 1);

            if (trim)
            {
                ret = ret.Trim();
            }

            return ret;
        }
    }

    public class info_stride
    {
    

        public abstract class stride_record
        {
            //1-3	Record code
            //4-5	Not used
            //6-73	Data
            //74-75	Not used
            //75-79	Four letter PDB code (if available)

            public string RecordCode;
            public string NotUsed1;
            public string Data;
            public string NotUsed2;
            public string PdbCode;

            protected stride_record(string column_format_line)
            {
                if (string.IsNullOrWhiteSpace(column_format_line)) return;

                RecordCode = column_format_line.column_value(1, 3);
                NotUsed1 = column_format_line.column_value(4, 5);
                Data = column_format_line.column_value(6, 73);
                NotUsed2 = column_format_line.column_value(74, 75);
                PdbCode = column_format_line.column_value(75, 79);
            }
        }

        public class Stride_Remark : stride_record
        {
            //REM Remarks and blank   lines

            public Stride_Remark(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_Header : stride_record
        {
            //HDR Header.  Protein name, date of file creation and PDB code

            public Stride_Header(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_Compound : stride_record
        {
            //CMP Compound.Full name of the molecule and identifying information

            public Stride_Compound(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_Receptor : stride_record
        {
            //SRC Species, organ, tissue, and mutant from which the molecule has been obtained

            public Stride_Receptor(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_StructureAuthors : stride_record
        {
            //AUT Names of the structure authors

            public Stride_StructureAuthors(string columnFormatLine) : base(columnFormatLine)
            {

            }

        }

        public class Stride_Chain : stride_record
        {
            //CHN File name and PDB chain identifier*).
            // Format: File name beginning from position 6 followed by one space and one-letter chain    identifier

            // 0        1         2         3         4         5         6         7
            // 1234567890123456789012345678901234567890123456789012345678901234567890123456789
            // CHN  /webclu/data/stride/pdb/pdb1a12.ent C                                 1A12

            public string FileName;
            public string OneLetterChainIdentifier;

            public Stride_Chain(string columnFormatLine) : base(columnFormatLine)
            {
                FileName = columnFormatLine.column_value(6, 75);
                OneLetterChainIdentifier = FileName.Substring(FileName.LastIndexOf(' ') + 1);
                FileName = FileName.Substring(0, FileName.LastIndexOf(' '));
            }
        }

        public class Stride_AminoAcidSequence : stride_record
        {
            //SEQ Amino acid sequence
            //
            //Format:
            //  6-9  First residue PDB number
            // 11-60 Sequence
            // 62-65 Last residue PDB number

            public string FirstResiduePdbNumber;
            public string Sequence;
            public string LastResiduePdbNumber;

            public Stride_AminoAcidSequence(string columnFormatLine) : base(columnFormatLine)
            {
                FirstResiduePdbNumber = columnFormatLine.column_value(6, 9);
                Sequence = columnFormatLine.column_value(11, 60);
                LastResiduePdbNumber = columnFormatLine.column_value(62, 65);
            }
        }

        public class Stride_SecondaryStructure : stride_record
        {

            //STR Secondary   structure summary
            //
            //Format:
            // 11-60 Secondary structure assignment**

            public string SecondaryStructureAssignment;

            public Stride_SecondaryStructure(string columnFormatLine) : base(columnFormatLine)
            {
                SecondaryStructureAssignment = columnFormatLine.column_value(11, 60);
            }
        }

        public class Stride_LocationOfSecondaryStructureElements : stride_record
        {
            //LOC Location of secondary structure elements
            //
            //Format:
            //  6-17 Element name
            // 19-21 First residue name
            // 22-27 First residue PDB number (spec says 32-26 !)
            // 28-28 First residue chain identifier
            // 36-38 Last residue name
            // 42-45 Last residue PDB number
            // 47-47 Last residue chain identifier

            // 0        1         2         3         4         5         6         7
            // 1234567890123456789012345678901234567890123456789012345678901234567890123456789
            // LOC  GammaInv     LEU   117 C      GLU    119 C                            1A12
            public string ElementName;
            public string FirstResidueName;
            public string FirstResiduePdbNumber;
            public string FirstResidueChainIdentifier;
            public string LastResidueName;
            public string LastResiduePdbNumber;
            public string LastResidueChainIdentifier;

            public Stride_LocationOfSecondaryStructureElements(string columnFormatLine) : base(columnFormatLine)
            {
                ElementName = columnFormatLine.column_value(6, 17);
                FirstResidueName = columnFormatLine.column_value(19, 21);
                FirstResiduePdbNumber = columnFormatLine.column_value(22, 27);
                FirstResidueChainIdentifier = columnFormatLine.column_value(28, 28);
                LastResidueName = columnFormatLine.column_value(36, 38);
                LastResiduePdbNumber = columnFormatLine.column_value(42, 45);
                LastResidueChainIdentifier = columnFormatLine.column_value(47, 47);
            }
        }

        public class Stride_DetailedSecondaryStructureAssignments : stride_record
        {
            //ASG Detailed secondary structure assignment
            //
            //Format:
            //  6-8  Residue name
            // 10-10 Protein chain identifier
            // 12-15 PDB residue number
            // 17-20 Ordinal residue number
            // 25-25 One letter secondary structure code**)
            // 27-39 Full secondary structure name
            // 43-49 Phi angle
            // 53-59 Psi angle
            // 65-69 Residue solvent accessible area

            public const char stride_default_secondary_structure = 'C';

            public string ResidueName;
            public string ProteinChainIdentifier;
            public string PdbResidueNumber;
            public string OrdinalResidueNumber;
            public string OneLetterSecondaryStructureCode;
            public string FullSecondaryStructureName;
            public string PhiAngle;
            public string PsiAngle;
            public string ResidueSolventAccessibleArea;

            public Stride_DetailedSecondaryStructureAssignments(string columnFormatLine) : base(columnFormatLine)
            {
                ResidueName = columnFormatLine.column_value(6, 8);
                ProteinChainIdentifier = columnFormatLine.column_value(10, 10);
                PdbResidueNumber = columnFormatLine.column_value(12, 15);
                OrdinalResidueNumber = columnFormatLine.column_value(17, 20);
                OneLetterSecondaryStructureCode = columnFormatLine.column_value(25, 25);
                FullSecondaryStructureName = columnFormatLine.column_value(27, 39);
                PhiAngle = columnFormatLine.column_value(43, 49);
                PsiAngle = columnFormatLine.column_value(53, 59);
                ResidueSolventAccessibleArea = columnFormatLine.column_value(65, 69);

                if (string.IsNullOrWhiteSpace(OneLetterSecondaryStructureCode)) OneLetterSecondaryStructureCode = stride_default_secondary_structure.ToString();
            }

        }

        public class Stride_LigandDonorReside : stride_record
        {
            // DNR LigandDonor residue
            //
            //Format:
            // 6-8  LigandDonor residue name
            //10-10 Protein chain identifier
            //12-15 PDB residue number
            //17-20 Ordinal residue number
            //26-28 Acceptor residue name
            //30-30 Protein chain identifier
            //32-35 PDB residue number
            //37-40 Ordinal residue number
            //42-45 N__0 distance
            //47-52 N__O_C angle
            //54-59 O__N_C angle
            //61-66 Angle between the planes of ligand_donor complex and O__N_C
            //68-73 Angle between the planes of acceptor complex and N__O_C

            public string LigandDonorResidueName;
            public string ProteinChainIdentifier1;
            public string PdbResidueNumber1;
            public string OrdinalResidueNumber1;
            public string AcceptorResidueName;
            public string ProteinChainIdentifier2;
            public string PdbResidueNumber2;
            public string OrdinalResidueNumber2;
            public string N__0Distance;
            public string N__o_cAngle;
            public string O__n_cAngle;
            public string AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c;
            public string AngleBetweenThePlanesOfAcceptorComplexAndN__o_c;

            public Stride_LigandDonorReside(string columnFormatLine) : base(columnFormatLine)
            {
                LigandDonorResidueName = columnFormatLine.column_value(6, 8);
                ProteinChainIdentifier1 = columnFormatLine.column_value(10, 10);
                PdbResidueNumber1 = columnFormatLine.column_value(12, 15);
                OrdinalResidueNumber1 = columnFormatLine.column_value(17, 20);
                AcceptorResidueName = columnFormatLine.column_value(26, 28);
                ProteinChainIdentifier2 = columnFormatLine.column_value(30, 30);
                PdbResidueNumber2 = columnFormatLine.column_value(32, 35);
                OrdinalResidueNumber2 = columnFormatLine.column_value(37, 40);
                N__0Distance = columnFormatLine.column_value(42, 45);
                N__o_cAngle = columnFormatLine.column_value(47, 52);
                O__n_cAngle = columnFormatLine.column_value(54, 59);
                AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c = columnFormatLine.column_value(61, 66);
                AngleBetweenThePlanesOfAcceptorComplexAndN__o_c = columnFormatLine.column_value(68, 73);
            }
        }

        public class Stride_AcceptorResidue : stride_record
        {
            // ACC Acceptor residue
            //
            //Format:
            // 6-8  Acceptor residue name
            //10-10 Protein chain identifier
            //12-15 PDB residue number
            //17-20 Ordinal residue number
            //26-28 LigandDonor residue name
            //30-30 Protein chain identifier
            //32-35 PDB residue number
            //37-40 Ordinal residue number
            //42-45 N..0 distance
            //47-52 N..O=C angle
            //54-59 O..N-C angle
            //61-66 Angle between the planes of ligand_donor complex and O..N-C
            //68-73 angle between the planes of acceptor complex and N..O=C

            public string AcceptorResidueName;
            public string ProteinChainIdentifier1;
            public string PdbResidueNumber1;
            public string OrdinalResidueNumber1;
            public string LigandDonorResidueName;
            public string ProteinChainIdentifier2;
            public string PdbResidueNumber2;
            public string OrdinalResidueNumber2;
            public string N__0Distance;
            public string N__o_cAngle;
            public string O__n_cAngle;
            public string AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c;
            public string AngleBetweenThePlanesOfAcceptorComplexAndN__o_c;

            public Stride_AcceptorResidue(string columnFormatLine) : base(columnFormatLine)
            {
                AcceptorResidueName = columnFormatLine.column_value(6, 8);
                ProteinChainIdentifier1 = columnFormatLine.column_value(10, 10);
                PdbResidueNumber1 = columnFormatLine.column_value(12, 15);
                OrdinalResidueNumber1 = columnFormatLine.column_value(17, 20);
                LigandDonorResidueName = columnFormatLine.column_value(26, 28);
                ProteinChainIdentifier2 = columnFormatLine.column_value(30, 30);
                PdbResidueNumber2 = columnFormatLine.column_value(32, 35);
                OrdinalResidueNumber2 = columnFormatLine.column_value(37, 40);
                N__0Distance = columnFormatLine.column_value(42, 45);
                N__o_cAngle = columnFormatLine.column_value(47, 52);
                O__n_cAngle = columnFormatLine.column_value(54, 59);
                AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c = columnFormatLine.column_value(61, 66);
                AngleBetweenThePlanesOfAcceptorComplexAndN__o_c = columnFormatLine.column_value(68, 73);
            }
        }


        public static List<stride_record> Load(string strideFilename)
        {
            var result = new List<stride_record>();

            var lines = io.ReadAllLines(strideFilename);

            foreach (var columnFormatLine in lines)
            {
                if (string.IsNullOrWhiteSpace(columnFormatLine) || columnFormatLine.Length <= 3) continue;

                var recordCode = columnFormatLine.Substring(0, 3);
                stride_record record = null;

                switch (recordCode)
                {
                    case "REM":
                        record = new Stride_Remark(columnFormatLine); break;
                    case "HDR":
                        record = new Stride_Header(columnFormatLine); break;
                    case "CMP":
                        record = new Stride_Compound(columnFormatLine); break;
                    case "SRC":
                        record = new Stride_Receptor(columnFormatLine); break;
                    case "AUT":
                        record = new Stride_StructureAuthors(columnFormatLine); break;
                    case "CHN":
                        record = new Stride_Chain(columnFormatLine); break;
                    case "SEQ":
                        record = new Stride_AminoAcidSequence(columnFormatLine); break;
                    case "STR":
                        record = new Stride_SecondaryStructure(columnFormatLine); break;
                    case "LOC":
                        record = new Stride_LocationOfSecondaryStructureElements(columnFormatLine); break;
                    case "ASG":
                        record = new Stride_DetailedSecondaryStructureAssignments(columnFormatLine); break;
                    case "DNR":
                        record = new Stride_LigandDonorReside(columnFormatLine); break;
                    case "ACC":
                        record = new Stride_AcceptorResidue(columnFormatLine); break;
                    default:
                        //record = new Stride_Record(columnFormatLine);
                        break;
                }

                if (record != null) result.Add(record);
            }

            return result;
        }

    }
}
