/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.maventest.eindopdrachtcourse5_1;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import javax.swing.JFileChooser;

/** The main GUI for simple nucleotide fasta analysis
 *
 * @author Thijs Weenink
 * @version 1.0
 */
public class EindopdrachtGUI extends javax.swing.JFrame {

    /**
     * Creates new form EindopdrachtGUI
     */
    public EindopdrachtGUI() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        fileChooser = new javax.swing.JFileChooser();
        jScrollPane1 = new javax.swing.JScrollPane();
        inputTextArea = new javax.swing.JTextArea();
        jLabel1 = new javax.swing.JLabel();
        errorLabel = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        headerField = new javax.swing.JTextField();
        jMenuBar1 = new javax.swing.JMenuBar();
        jMenu1 = new javax.swing.JMenu();
        openFileMenu = new javax.swing.JMenuItem();
        analyseMenu = new javax.swing.JMenu();
        nuclPercentagesMenu = new javax.swing.JMenuItem();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Eindopdracht Course 5");

        inputTextArea.setColumns(20);
        inputTextArea.setRows(5);
        jScrollPane1.setViewportView(inputTextArea);

        jLabel1.setFont(new java.awt.Font("Tahoma", 0, 18)); // NOI18N
        jLabel1.setText("Input:");

        errorLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);

        jLabel2.setFont(new java.awt.Font("Tahoma", 0, 18)); // NOI18N
        jLabel2.setText("Header:");

        jMenu1.setText("File");

        openFileMenu.setText("Open...");
        openFileMenu.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                openFileMenuActionPerformed(evt);
            }
        });
        jMenu1.add(openFileMenu);

        jMenuBar1.add(jMenu1);

        analyseMenu.setText("Analyse");

        nuclPercentagesMenu.setText("Calc nulc. percentages");
        nuclPercentagesMenu.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                nuclPercentagesMenuActionPerformed(evt);
            }
        });
        analyseMenu.add(nuclPercentagesMenu);

        jMenuBar1.add(analyseMenu);

        setJMenuBar(jMenuBar1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(headerField)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(errorLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 450, Short.MAX_VALUE)
                            .addComponent(jLabel1, javax.swing.GroupLayout.Alignment.LEADING))
                        .addComponent(jLabel2)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(headerField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(errorLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 70, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void openFileMenuActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openFileMenuActionPerformed
        //Opening of the file via fileChooser
        int returnVal = fileChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            try {
                // What to do with the file, e.g. display it in a TextArea
                inputTextArea.read( new FileReader( file.getAbsolutePath() ), null );
                
                //Displaying the header in the textField
                String text = inputTextArea.getText();
                for(String line:text.split("\n")){
                    if(line.startsWith(">")){
                        headerField.setText(line);
                    }
                } 
            } catch (IOException ex) {
                System.out.println("problem accessing file"+file.getAbsolutePath());
            }
        } else {
            System.out.println("File access cancelled by user.");
        }
               
    }//GEN-LAST:event_openFileMenuActionPerformed

    private void nuclPercentagesMenuActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_nuclPercentagesMenuActionPerformed
        /*
        Simple logic for displaying a PieChart for the nucleotide percentages.
        All the methods used are written in the classes: AminoUtils.java and PieChart.java.
        AminoUtils is a collection of methods. Only 2 are used in total. (Bottom 2 methods)
        
        Only if the try{}catch{} gets passed without errors, the PieChart will be made.
        */
        
        String fixedFasta;
        String fastaHeader = null;
        int fastaLength = 0;
        boolean noError = true;
        try {    
            fixedFasta = AminoUtils.removeBreaklines(inputTextArea.getText());
            nuclPercentagesArray = AminoUtils.getNuclPercentages(fixedFasta.toUpperCase());
            fastaHeader = fixedFasta.split("\n")[0];
            fastaLength = (fixedFasta.split("\n")[1]).length();
        } catch (IncorrectFileFormatError ex) {
            noError = false;
            errorLabel.setText(ex.toString());
            System.out.println(ex.toString());
        } catch (UnknownNucleotideError ex) {
            noError = false;
            errorLabel.setText(ex.toString());
            System.out.println(ex.toString());
        } catch (ArrayIndexOutOfBoundsException ex){
            noError = false;
            errorLabel.setText(ex.toString());
            System.out.println(ex.toString());
        } catch (Exception ex){
            noError = false;
            errorLabel.setText(ex.toString());
            System.out.println(ex.toString());
        }
        
        if(noError){
            String[] nucls = {"A", "T", "C", "G"};
            PieChart.showChart(nucls, nuclPercentagesArray, "Piechart of "+fastaHeader, fastaLength);
        }
        
    }//GEN-LAST:event_nuclPercentagesMenuActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(EindopdrachtGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(EindopdrachtGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(EindopdrachtGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(EindopdrachtGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new EindopdrachtGUI().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenu analyseMenu;
    private javax.swing.JLabel errorLabel;
    private javax.swing.JFileChooser fileChooser;
    private javax.swing.JTextField headerField;
    private javax.swing.JTextArea inputTextArea;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JMenu jMenu1;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JMenuItem nuclPercentagesMenu;
    private javax.swing.JMenuItem openFileMenu;
    // End of variables declaration//GEN-END:variables
    private double[] nuclPercentagesArray;

}
