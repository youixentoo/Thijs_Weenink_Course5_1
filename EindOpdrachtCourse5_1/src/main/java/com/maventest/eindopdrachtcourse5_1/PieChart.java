package com.maventest.eindopdrachtcourse5_1;


import java.awt.BorderLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.general.DefaultPieDataset;
import org.jfree.data.general.PieDataset;

/** Pop-up frame with a PieChart from JFreeChart
 *
 * @author Thijs Weenink
 * @version 1.0
 */
public class PieChart extends JFrame {
    
    private PieChart() {
        /*
        All the variables are declared and assigned outside this method
        */
        PieDataset dataset = createDataset();
        JPanel chartPanel = createChartPanel(dataset);
        super.setTitle(chartTitle);
        super.add(chartPanel, BorderLayout.CENTER);
 
        super.setSize(xWidth, yWidth);
        super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        super.setLocationRelativeTo(null);
    }
    
    private JPanel createChartPanel(PieDataset dataset){
        /*
        Creates the piechart which will be displayed on the frame
        The title gets used to display the length of the sequence
        */
        JFreeChart chart = ChartFactory.createPieChart("Total length of sequence is: "+totalLength, dataset);
        
        return new ChartPanel(chart);
    }
    
    private PieDataset createDataset(){
        /*
        Creates the dataset which will be displayed in the piechart
        The arrayIndex pairs the nucleotide with its percentage
        */
        
        DefaultPieDataset dataset = new DefaultPieDataset();
        int arrayIndex = 0;
              
        for (double nPerc:nuclPercs){
            dataset.setValue(nucleotides[arrayIndex], nPerc);
            arrayIndex++;
        }

        return dataset;
    }
    
    /** Creates a pop-up frame with a Nucleotide piechart.
     * <p> Hover of the different pieces to view the percentages.
     *
     * @param nucleotides A String[] with the used nucleotides, in the same order as nuclPercentages
     * @param nuclPercentages A double[] with the percentages of the nucleotides, in the same order as nucleotides
     * @param chartTitle The title of the graph
     * @param sequenceLength The length of the sequence
     */
    public static void showChart(String[] nucleotides, double[] nuclPercentages, String chartTitle, int sequenceLength){
        // Variable assignment
        nuclPercs = nuclPercentages;
        PieChart.nucleotides = nucleotides;
        PieChart.chartTitle = chartTitle;
        totalLength = sequenceLength;
        
        SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
            new PieChart().setVisible(true);
        }
        });
        
    }
    // Variable declaration
    private static String chartTitle;
    private static int xWidth = 640;
    private static int yWidth = 480;
    private static double[] nuclPercs;
    private static String[] nucleotides;
    private static int totalLength;
}