/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.maventest.eindopdrachtcourse5_1;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * A collection of DNA, RNA and protein utilities.<br>
 * Translations, transcriptions, nucleotide analysis, protein analysis, protein codes conversion, linebreak removal.
 * 
 * @author Thijs Weenink
 * @version 1.2
 */
public class AminoUtils {
    
    /**
     *
     * @param seq The DNA sequence as a String
     * @return The peptide sequence in one codes.
     * <p>Returns "No full protein available" if the DNA sequence does not give a full protein.
     * <br>This gives no exception so it can be used easily with a multiple fasta formatted file.
     * @throws UnknownNucleotideError If there is letter other then ATGC in the sequence.
     * @throws StringIndexOutOfBoundsException If there is no protein to be gotten
     * <p>The StringIndexOutOfBoundsException will most likely only be thrown when used in a for loop when reading a multiple fasta file.
     */
    public static String getTranslation(String seq) throws UnknownNucleotideError{
            String protein = "";
            Integer start_Codon;

            java.util.Map table = new HashMap(); //Codon table declaration
            // <editor-fold defaultstate="collapsed" desc="Table contents">
            table.put ("TTT", "F");
            table.put ("TTC", "F");
            table.put ("TTA", "L");
            table.put ("TTG", "L");
            table.put ("TCT", "S");
            table.put ("TCC", "S");
            table.put ("TCA", "S");
            table.put ("TCG", "S");
            table.put ("TAT", "Y");
            table.put ("TAC", "Y");
              //            TAA end
              //            TAG end
            table.put ("TGT", "C");
            table.put ("TGC", "C");
              //            TGA end
            table.put ("TGG", "W");
            table.put ("CTT", "L");
            table.put ("CTC", "L");
            table.put ("CTA", "L");
            table.put ("CTG", "L");
            table.put ("CCT", "P");
            table.put ("CCC", "P");
            table.put ("CCA", "P");
            table.put ("CCG", "P");
            table.put ("CAT", "H");
            table.put ("CAC", "H");
            table.put ("CAA", "Q");
            table.put ("CAG", "Q");
            table.put ("CGT", "R");
            table.put ("CGC", "R");
            table.put ("CGA", "R");
            table.put ("CGG", "R");
            table.put ("ATT", "I");
            table.put ("ATC", "I");
            table.put ("ATA", "I");
            table.put ("ATG", "M");
            table.put ("ACT", "T");
            table.put ("ACC", "T");
            table.put ("ACA", "T");
            table.put ("ACG", "T");
            table.put ("AAT", "N");
            table.put ("AAC", "N");
            table.put ("AAA", "K");
            table.put ("AAG", "K");
            table.put ("AGT", "S");
            table.put ("AGC", "S");
            table.put ("AGA", "R");
            table.put ("AGG", "R");
            table.put ("GTT", "V");
            table.put ("GTC", "V");
            table.put ("GTA", "V");
            table.put ("GTG", "V");
            table.put ("GCT", "A");
            table.put ("GCC", "A");
            table.put ("GCA", "A");
            table.put ("GCG", "A");
            table.put ("GAT", "D");
            table.put ("GAC", "D");
            table.put ("GAA", "E");
            table.put ("GAG", "E");
            table.put ("GGT", "G");
            table.put ("GGC", "G");
            table.put ("GGA", "G");
            table.put ("GGG", "G");
            // Source: http://gaggle.systemsbiology.net/svn/gaggle/gaggle/branches/original/org/systemsbiology/gaggle/geese/sequence/CodonTable.java
            // </editor-fold>

            if(seq.matches(".*[ATCG].*")){
                start_Codon = seq.indexOf("ATG");          
                String string_aug = seq.substring(start_Codon,seq.length());
                String[] splitted_codons = string_aug.split("(?<=\\G...)");

                for (String codon:splitted_codons){ //Translation loop
                    if(codon.equals("TAA")){
                        break;
                    }else if(codon.equals("TAG")){
                        break;
                    }else if(codon.equals("TGA")){
                        break;
                    }else if(codon.length() < 3){//The way the string is splitted means there can be "codons" with a length less than 3, if this is true, there is no (full) protein.
                        protein = "No full protein available";
                        break;
                    }else{
                        protein += (String)table.get(codon);
                    }
                }

                return protein;
            }else{
                throw new UnknownNucleotideError("UnknownNucleotideError Occurred:\nThere is an unknown nucleotide in this sequence.");
            }
            
    }
    
    /**
     *
     * @param seq The one letter peptide sequence as a String
     * @param fullName A boolean to indicate whether or not the fullnames will be returned. <br>true = fullnames<br>false = three letter code.
     * @return Either the three letter code or the fullnames (String).
     * @throws UnknownAminoacidException If there is an unknown aminoacid in the sequence.
     */
    public static String getsequenceThree(String seq, boolean fullName) throws UnknownAminoacidException{
        int error_index;
        String threeletterCode = "";
        java.util.Map table = new HashMap(); //Codon table declaration
        // <editor-fold defaultstate="collapsed" desc="3 letter Aminoacids">
        table.put ("F", "Phe"); //1
        table.put ("L", "Leu"); //2
        table.put ("A", "Ala"); //3
        table.put ("R", "Arg"); //4
        table.put ("N", "Asn"); //5
        table.put ("D", "Asp"); //6
        table.put ("C", "Cys"); //7
        table.put ("Q", "Gln"); //8
        table.put ("E", "Glu"); //9
        table.put ("G", "Gly"); //10
        table.put ("H", "His"); //11
        table.put ("I", "Ile"); //12
        table.put ("K", "Lys"); //13
        table.put ("M", "Met"); //14
        table.put ("P", "Pro"); //15
        table.put ("S", "Ser"); //16 
        table.put ("T", "Thr"); //17
        table.put ("W", "Trp"); //18
        table.put ("Y", "Tyr"); //19
        table.put ("V", "Val"); //20
        // </editor-fold>
        
        // <editor-fold defaultstate="collapsed" desc="Fullnames Aminoacids">
        java.util.Map table2 = new HashMap();
        table2.put ("A", "Alenine");
        table2.put ("R", "Arginine");
        table2.put ("N", "Asparagine");
        table2.put ("C", "Cysteine");
        table2.put ("D", "Aspartate");
        table2.put ("Q", "Glutamine");
        table2.put ("E", "Glutamate");
        table2.put ("G", "Glycine");
        table2.put ("H", "Histidine");
        table2.put ("I", "Isoleucine");
        table2.put ("L", "Leucine");
        table2.put ("K", "Lysine");
        table2.put ("M", "Methionine");
        table2.put ("F", "Phenylalanine");
        table2.put ("P", "Proline");
        table2.put ("S", "Serine");
        table2.put ("T", "Threonine");
        table2.put ("W", "Thryptophan");
        table2.put ("V", "Valine");
        table2.put ("Y", "Tyrosine");
        // </editor-fold>
            
        String[] sequence = seq.split("(?<=\\G.)");
        for (String letter:sequence){
            if(fullName){
                if((String)table2.get(letter.toUpperCase()) != null){
                    threeletterCode += (String)table2.get(letter.toUpperCase());
                }else{
                    error_index = ((threeletterCode.split("(?=\\p{Upper})")).length)+1;
                    System.out.println(error_index);
                    threeletterCode = null;
                    throw new UnknownAminoacidException("UnknownAminoacidException Occurred:\nThere is an unknown 1 letter aminoacid in this sequence, at location: "+error_index);
                }
            }else{
                if((String)table.get(letter.toUpperCase()) != null){
                    threeletterCode += (String)table.get(letter.toUpperCase());
                }else{
                    error_index = (threeletterCode.length()/3)+1;
                    System.out.println("Error index: "+error_index);
                    threeletterCode = null;
                    throw new UnknownAminoacidException("UnknownAminoacidException Occurred:\nThere is an unknown 1 letter aminoacid in this sequence, at location: "+error_index);
                }
            }
        }
        return threeletterCode;
    }
    
    /**
     *
     * @param seq The three letter peptide sequence as a String
     * @return The one letter peptide sequence (String)
     * @throws UnknownAminoacidException If there is an unknown aminoacid in the sequence.
     */
    public static String getsequenceOne(String seq) throws UnknownAminoacidException{
        int error_index;
        String oneletterCode = "";
         java.util.Map table = new HashMap(); //Codon table declaration
            // <editor-fold defaultstate="collapsed" desc="Table contents">
            table.put ("Phe", "F");
            table.put ("Leu", "L");
            table.put ("Ala", "A");
            table.put ("Arg", "R");
            table.put ("Asn", "N");
            table.put ("Asp", "D");
            table.put ("Cys", "C");
            table.put ("Gln", "Q");
            table.put ("Glu", "E");
            table.put ("Gly", "G");
            table.put ("His", "H");
            table.put ("Ile", "I");
            table.put ("Lys", "K");
            table.put ("Met", "M");
            table.put ("Pro", "P");
            table.put ("Ser", "S");
            table.put ("Thr", "T");
            table.put ("Trp", "W");
            table.put ("Tyr", "Y");
            table.put ("Val", "V");
            // </editor-fold>
        
        
        
        for (String Abr:seq.split("(?<=\\G...)")){
            if((String)table.get(Abr) != null){
                oneletterCode += (String)table.get(Abr);
            }else{
                error_index = oneletterCode.length()+1;
                oneletterCode = null;
                throw new UnknownAminoacidException("UnknownAminoacidException Occurred:\nThere is an unknown 3 letter aminoacid in this sequence, at location: "+error_index);
            }   
        }
        return oneletterCode;
    }
    
    /**
     *
     * @param seq The one letter peptide sequence as a String 
     * @return A double[] which contains the two polarity percentages
     * <p>[0] = The percentage of polar aminoacids in the sequence.
     * <p>[1] = The percenatge of apolar aminoacids in the sequence.
     * @throws UnknownAminoacidException If there is an unknown aminoacid in the sequence.
     */
    public static double[] getPolarity(String seq) throws UnknownAminoacidException{
        java.util.HashMap polarity_table = new HashMap();
        // <editor-fold defaultstate="collapsed" desc="Polarity Table">
        polarity_table.put("A", 0.0);
        polarity_table.put("R", 1.0);
        polarity_table.put("N", 1.0);
        polarity_table.put("D", 1.0);
        polarity_table.put("C", 1.0);
        polarity_table.put("F", 0.0);
        polarity_table.put("Q", 1.0);
        polarity_table.put("E", 1.0);
        polarity_table.put("G", 1.0);
        polarity_table.put("H", 1.0);
        polarity_table.put("I", 0.0);
        polarity_table.put("L", 0.0);
        polarity_table.put("K", 1.0);
        polarity_table.put("M", 0.0);
        polarity_table.put("P", 0.0);
        polarity_table.put("S", 1.0);
        polarity_table.put("T", 1.0);
        polarity_table.put("W", 0.0);
        polarity_table.put("Y", 1.0);
        polarity_table.put("V", 0.0); // 20 aminoacids with if they are polar (1) or not (0)
        // </editor-fold>
        
        double polarity_number = 0.0;
        String[] aminoArray = seq.split("(?<=\\G.)");
        double polar_perc, apolar_perc;
        
        
        
        for(String Acid:aminoArray){
            if(!Double.isNaN((double)polarity_table.get(Acid))){
                polarity_number += (double)polarity_table.get(Acid);
            }else{
                throw new UnknownAminoacidException("UnknownAminoacidException Occurred:\nThere is an unknown aminoacid in this sequence.");
            }
        }
        
        polar_perc = (polarity_number / aminoArray.length)*100;
        apolar_perc = 100 - polar_perc;
        double[] polarity_array = {polar_perc, apolar_perc};
        
        //System.out.println(polarity_number+" "+aminoArray.length+" "+polar_perc+" "+apolar_perc);

        return polarity_array;   
    }
    
    /**
     *
     * @param seq The one letter peptide sequence (String)
     * @return A double[] which contains the 4 values.
     * <p>[0] = Charge number (double), a positive number indicates that the protein contains more positive aminoacids then negative, negative the other way around and 0.0 that the number of negative and positive aminoacids are the same.
     * <p>[1] = positive (double), the amount of positive aminoacids in the sequence.
     * <p>[2] = neutral (double), the amount of neutral aminoacids in the sequence.
     * <p>[3] = negative (double), the amount of negative aminoacids in the sequence.
     * @throws UnknownAminoacidException If there is an unknown aminoacid in the sequence.
     */
    public static double[] getCharges(String seq) throws UnknownAminoacidException{
        java.util.HashMap charge_table = new HashMap();
        //<editor-fold defaultstate="collapsed" desc="Charge table">
        charge_table.put("A", 0.0);
        charge_table.put("R", 1.0);
        charge_table.put("N", 0.0);
        charge_table.put("D", -1.0);
        charge_table.put("C", 0.0);
        charge_table.put("F", 0.0);
        charge_table.put("Q", 0.0);
        charge_table.put("E", -1.0);
        charge_table.put("G", 0.0);
        charge_table.put("H", 1.0);
        charge_table.put("I", 0.0);
        charge_table.put("L", 0.0);
        charge_table.put("K", 1.0);
        charge_table.put("M", 0.0);
        charge_table.put("P", 0.0);
        charge_table.put("S", 0.0);
        charge_table.put("T", 0.0);
        charge_table.put("W", 0.0);
        charge_table.put("Y", 0.0);
        charge_table.put("V", 0.0); // 20 aminoacids with if they are positive (1), neutral(0) and negative (-1)
        //</editor-fold>
    
        double charge_number = 0.0;
        String[] aminoArray = seq.split("(?<=\\G.)");
        double positive = 0, neutral = 0, negative = 0;
        
        for(String Acid:aminoArray){
            if(!Double.isNaN((double)charge_table.get(Acid))){
                charge_number += (double)charge_table.get(Acid);
                if((double)charge_table.get(Acid) == -1.0){
                    negative++;
                }else if((double)charge_table.get(Acid) == 1.0){
                    positive++;
                }else{
                    neutral++;
                }
            }else{
                throw new UnknownAminoacidException("UnknownAminoacidException Occurred:\nThere is an unknown aminoacid in this sequence.");
            }
        }
        
        double[] charge_Array = {charge_number, positive, neutral, negative};

        return charge_Array;
    }
    
    /**
     *
     * @param seq The one letter peptide sequence. (String)
     * @return List of doubles containing the Hydropathy values.
     * @throws NotCodingPeptideException If the sequence does not start with "M".
     * @throws UnknownAminoacidException If the sequence contains an unknown aminoacid.
     */
    public static List<Double> getHydropathy(String seq) throws NotCodingPeptideException, UnknownAminoacidException{
        
        java.util.HashMap hydro_table = new HashMap();
        //<editor-fold defaultstate="collapsed" desc="Hydropaty table">
        hydro_table.put("I",-0.528);
        hydro_table.put("L",-0.342);
        hydro_table.put("F",-0.370);
        hydro_table.put("V",-0.308);
        hydro_table.put("M",-0.324);
        hydro_table.put("P",-0.322);
        hydro_table.put("W",-0.270);
        hydro_table.put("H",2.029);
        hydro_table.put("T",0.853);
        hydro_table.put("E",3.173);
        hydro_table.put("Q",2.176);
        hydro_table.put("C",0.081);
        hydro_table.put("Y",1.677);
        hydro_table.put("A",-0.495);
        hydro_table.put("S",0.936);
        hydro_table.put("N",2.354);
        hydro_table.put("D",9.573);
        hydro_table.put("R",4.383);
        hydro_table.put("G",0.386);
        hydro_table.put("K",2.101);
        //</editor-fold>
        
        if(!seq.startsWith("M")){
            throw new NotCodingPeptideException("NotPeptideException Occurred:\nThis sequence doesnt start with Methione (M)");
        }
        
        List<Double> HydroArrayList = new ArrayList<>(seq.length());
        String[] aminoArray = seq.split("(?<=\\G.)");
        
        for(String Acid:aminoArray){
            if(!Double.isNaN((double)hydro_table.get(Acid))){
                HydroArrayList.add((double)hydro_table.get(Acid));
            }else{
                throw new UnknownAminoacidException("UnknownAminoacidException Occurred:\nThere is an unknown aminoacid in this sequence.");
            }
        }        
        return HydroArrayList;
    }
    
    /**
     *
     * @param seq The one letter peptide sequence. (String)
     * @return List of "Pho", "Phi" and "Neu" for each aminoacid in the sequence.
     * @throws UnknownAminoacidException If the sequence contains an unknown aminoacid.
     */
    public static List<String> getHydrophobicities(String seq) throws UnknownAminoacidException{
        java.util.HashMap Phobicity_table = new HashMap(); // pH 7 values: http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html#hydro
        //<editor-fold defaultstate="collapsed" desc="Hydrophobicity table">
        Phobicity_table.put("L", "Pho");
        Phobicity_table.put("I", "Pho");
        Phobicity_table.put("F", "Pho");
        Phobicity_table.put("W", "Pho");
        Phobicity_table.put("V", "Pho"); // Hydrophobic
        Phobicity_table.put("M", "Pho");
        Phobicity_table.put("C", "Pho");
        Phobicity_table.put("Y", "Pho");
        Phobicity_table.put("A", "Pho");
        
        Phobicity_table.put("T", "Neu");
        Phobicity_table.put("H", "Neu");
        Phobicity_table.put("G", "Neu"); // neutral
        Phobicity_table.put("S", "Neu");
        Phobicity_table.put("Q", "Neu");
        
        Phobicity_table.put("R", "Phi");
        Phobicity_table.put("K", "Phi");
        Phobicity_table.put("N", "Phi"); // Hyrophilic
        Phobicity_table.put("E", "Phi");
        Phobicity_table.put("P", "Phi");
        Phobicity_table.put("D", "Phi");
        //</editor-fold>
        

        List<String> Phobicity_List = new ArrayList<>(seq.length() );
        String[] aminoArray = seq.split("(?<=\\G.)");
        
        for(String Acid:aminoArray){
            if((String)Phobicity_table.get(Acid) != null){
                Phobicity_List.add((String)Phobicity_table.get(Acid));
            }else{
                throw new UnknownAminoacidException("UnknownAminoacidException Occurred:\nThere is an unknown aminoacid in this sequence.");
            }
        }
        
        
        return Phobicity_List;
        
        
    }
    
    /**
     *
     * @param seq The DNA sequence as a String
     * @return The RNA sequence as a String
     * @throws UnknownNucleotideError If the sequence contains an unknown aminoacid.
     */
    public static String getRNA(String seq) throws UnknownNucleotideError{
        String RNA= null;
        if(seq.matches(".*[ATCG].*")){
            RNA = seq.replace("T", "U");   
        }else{
            throw new UnknownNucleotideError("UnknownNucleotideError Occurred:\nThere is an unknown nucleotide in this sequence.");
        }

        return RNA;
    }
    
    /**
     *
     * @param seq The DNA sequence as a String.
     * @return The protein sequence in one-letter codes (String)
     * @throws UnknownNucleotideError If there is a letter other than ATGC in the sequence
     * @throws NotAFullproteinError If the sequence is not dividable by 3.
     */
    public static String getprotein(String seq) throws UnknownNucleotideError, NotAFullproteinError{
            String protein = "";
            Integer start_Codon;
            seq = seq.toUpperCase();
            

            java.util.Map table = new HashMap(); //Codon table declaration
            // <editor-fold defaultstate="collapsed" desc="DNA -> protein table">
            table.put ("TTT", "F");
            table.put ("TTC", "F");
            table.put ("TTA", "L");
            table.put ("TTG", "L");
            table.put ("TCT", "S");
            table.put ("TCC", "S");
            table.put ("TCA", "S");
            table.put ("TCG", "S");
            table.put ("TAT", "Y");
            table.put ("TAC", "Y");
              //            TAA end
              //            TAG end
            table.put ("TGT", "C");
            table.put ("TGC", "C");
              //            TGA end
            table.put ("TGG", "W");
            table.put ("CTT", "L");
            table.put ("CTC", "L");
            table.put ("CTA", "L");
            table.put ("CTG", "L");
            table.put ("CCT", "P");
            table.put ("CCC", "P");
            table.put ("CCA", "P");
            table.put ("CCG", "P");
            table.put ("CAT", "H");
            table.put ("CAC", "H");
            table.put ("CAA", "Q");
            table.put ("CAG", "Q");
            table.put ("CGT", "R");
            table.put ("CGC", "R");
            table.put ("CGA", "R");
            table.put ("CGG", "R");
            table.put ("ATT", "I");
            table.put ("ATC", "I");
            table.put ("ATA", "I");
            table.put ("ATG", "M");
            table.put ("ACT", "T");
            table.put ("ACC", "T");
            table.put ("ACA", "T");
            table.put ("ACG", "T");
            table.put ("AAT", "N");
            table.put ("AAC", "N");
            table.put ("AAA", "K");
            table.put ("AAG", "K");
            table.put ("AGT", "S");
            table.put ("AGC", "S");
            table.put ("AGA", "R");
            table.put ("AGG", "R");
            table.put ("GTT", "V");
            table.put ("GTC", "V");
            table.put ("GTA", "V");
            table.put ("GTG", "V");
            table.put ("GCT", "A");
            table.put ("GCC", "A");
            table.put ("GCA", "A");
            table.put ("GCG", "A");
            table.put ("GAT", "D");
            table.put ("GAC", "D");
            table.put ("GAA", "E");
            table.put ("GAG", "E");
            table.put ("GGT", "G");
            table.put ("GGC", "G");
            table.put ("GGA", "G");
            table.put ("GGG", "G");
            // Source: http://gaggle.systemsbiology.net/svn/gaggle/gaggle/branches/original/org/systemsbiology/gaggle/geese/sequence/CodonTable.java
            // </editor-fold>

            if(seq.matches(".*[ATCG].*")){
                start_Codon = seq.indexOf("ATG"); 
                /*
                This can give a StringIndexOutOfBoundsException Exception, which means there is no protein to be gotten. 
                And is therefore useless, and the header and value gets skipped if thrown
                */         
                String string_aug = seq.substring(start_Codon,seq.length());
                String[] splitted_codons = string_aug.split("(?<=\\G...)");
                //System.out.println(header+"\n"+java.util.Arrays.toString(splitted_codons));


                for (String codon:splitted_codons){ //Translation loop
                    if(codon.equals("TAA")){
                        break;
                    }else if(codon.equals("TAG")){
                        break;
                    }else if(codon.equals("TGA")){
                        break;
                    }else if(codon.length() < 3){//The way the string is splitted means there can be "codons" with a length less than 3, if this is true, there is no (full) protein.
                        throw new NotAFullproteinError("<html><body>NotAFullproteinError Occurred:<br>This DNA sequence cannot be translated into a protein sequence.</body></html>");
                    }else{
                        protein += (String)table.get(codon);
                    }
                }
            }else{
                throw new UnknownNucleotideError("UnknownNucleotideError Occurred:\nThere is an unknown nucleotide in this sequence.");
            }
            
        return protein;
    }
    
    /**
     *
     * @param seq The DNA sequence as a String
     * @return The GC percentage as a double
     * @throws UnknownNucleotideError If there is a letter other than ATGC in the sequence
     */
    public static double getGCPercentage(String seq) throws UnknownNucleotideError{
        double gcAmount = 0;
        double gcPercentage;
        if(seq.matches(".*[ATCG].*")){
            for (String nucl:seq.split("(?<=\\G.)")){
                if("G".equals(nucl) || "C".equals(nucl)){
                    gcAmount++;
                }
            }  
        }else{
            throw new UnknownNucleotideError("UnknownNucleotideError Occurred:\nThere is an unknown nucleotide in this sequence.");
        }
        
        gcPercentage = (gcAmount/seq.length())*100;

        return gcPercentage;
    }
    
    /**
     *
     * @param seq The DNA sequence as a String
     * @return A double[] which contains the 4 values.
     * <br>[0] Contains the percentage of A in the sequence
     * <br>[1] Contains the percentage of T in the sequence
     * <br>[2] Contains the percentage of C in the sequence
     * <br>[3] Contains the percentage of G in the sequence
     * @throws UnknownNucleotideError If there is a letter other than ATGC in the sequence
     */
    public static double[] getNuclPercentages(String seq) throws UnknownNucleotideError{
        double nucl_A = 0; double nucl_T = 0; double nucl_C = 0; double nucl_G = 0;
        double perc_A, perc_T, perc_C, perc_G;
        String sequenceDNA = (seq.split("\n"))[1];
        if(sequenceDNA.matches(".*[ATCG].*")){
            
            for (String nucl:sequenceDNA.split("(?<=\\G.)")){
                if("A".equals(nucl)){
                    nucl_A++;
                }else if("T".equals(nucl)){
                    nucl_T++;
                }else if("C".equals(nucl)){
                    nucl_C++;
                }else if("G".equals(nucl)){
                    nucl_G++;
                }
            }
            
            perc_A = (nucl_A / sequenceDNA.length()) * 100;
            perc_T = (nucl_T / sequenceDNA.length()) * 100;
            perc_C = (nucl_C / sequenceDNA.length()) * 100;
            perc_G = (nucl_G / sequenceDNA.length()) * 100;

        }else{
            throw new UnknownNucleotideError("UnknownNucleotideError Occurred:\nThere is an unknown nucleotide in this sequence.");
        }
        
        double[] percentageArray = {perc_A, perc_T, perc_C, perc_G};
        
        return percentageArray;
    }
    
    /**
     *
     * @param fasta The fasta sequence with linebreaks.
     * @return A String with the linebreaks removed in the given sequence.
     * @throws IncorrectFileFormatError If the first character of the input does not match the greater than sign.
     */
    public static String removeBreaklines(String fasta) throws IncorrectFileFormatError{
        String[] lines = fasta.split("\\n");
        String fixed_Fasta = "";
        
        if(!lines[0].startsWith(">")){
            throw new IncorrectFileFormatError("IncorrectFileFormatError occurred:\nThe input does not start with a '>'");
        }
        
        for (String line: lines){
            if (line.startsWith(">")){
                fixed_Fasta += (line+"\n");
            }else if(!line.startsWith(">")){
                fixed_Fasta += (line);
            }
        } 
        return fixed_Fasta;
    }
    
    /*
    public String proteinsequence(String seq){
        String proteinsequence = "";
        int Index_M = seq.indexOf("M");
        
        
        
        return proteinsequence;
    }*/
        
}

/* Custom Exceptions */

class UnknownAminoacidException extends Exception{
    String str1;
    
    UnknownAminoacidException() {
        super();
    }

    public UnknownAminoacidException(String str2){
        super(str2);
        str1 = str2;
    }
    @Override
    public String toString() {
        System.out.println(str1);
        return (str1);
    }
}

class NotCodingPeptideException extends Exception{
    String str1;

    NotCodingPeptideException() {
        super();
    }
    
    public NotCodingPeptideException(String str2){
        super(str2);
        str1 = str2;
    }
    
    @Override
    public String toString(){
        System.out.println(str1);
        return str1;
    }  
}

class UnknownNucleotideError extends Exception{
    String str1;

    UnknownNucleotideError() {
        super();
    }

    public UnknownNucleotideError(String str2) {
        super(str2);
        str1 = str2;
    }
    
    @Override
    public String toString(){
        System.err.println(str1);
        return str1;
    }  
}

class NotAFullproteinError extends Exception{
    String str1;

    NotAFullproteinError() {
        super();
    }

    public NotAFullproteinError(String str2) {
        super(str2);
        str1 = str2;
    }
    
    @Override
    public String toString(){
        System.err.println(str1);
        return str1;  
    }   
}

class IncorrectFileFormatError extends Exception{
    String str1;

    IncorrectFileFormatError() {
        super();
    }

    public IncorrectFileFormatError(String str2) {
        super(str2);
        str1 = str2;
    }
    
    @Override
    public String toString(){
        System.err.println(str1);
        return str1;  
    }   
}


