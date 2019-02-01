import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class Toolbox {

    public static double[] averagedResults(double[][] inputData){

        int nReps = inputData.length;
        int nMeasurements = inputData[0].length;

        double[] averagedResults = new double[nMeasurements];

        //iterate over the counters, checking all the reps, then moving to next counter
        for(int c = 0; c < nMeasurements; c++){
            double runningTotal = 0.;
            for(int r = 0; r < nReps; r++){
                runningTotal += inputData[r][c];
            }
            averagedResults[c] = runningTotal/nReps;
        }
        return averagedResults;
    }

    public static double[][] averagedResults(double[][][] inputData){

        int nReps = inputData.length;
        int nTimes = inputData[0].length;
        int L = inputData[0][0].length;

        double[][] averagedResults = new double[nTimes][L];

        for(int t = 0; t < nTimes; t++){

            for(int l = 0; l < L; l++){

                double runningTotal = 0.;

                for(int r = 0; r < nReps; r++){
                    runningTotal += inputData[r][t][l];
                }
                averagedResults[t][l] = runningTotal/(double)nReps;
            }
        }
        return averagedResults;
    }


    public static void writeAveragedDistbsToFile(String filename, double[][] inputData){

        try {
            File file = new File(filename+".txt");
            if(!file.exists()) file.createNewFile();

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            int nTimes = inputData.length;
            int L = inputData[0].length;

            for(int l = 0; l < L; l++){
                String output = String.valueOf(l)+" ";
                for(int t = 0; t < nTimes; t++){
                    output += String.format("%.6f ", inputData[t][l]);
                }
                bw.write(output);
                bw.newLine();
            }
            bw.close();
        }catch (IOException e){}
    }




    public static void writeAveragedArrayToFile(String filename, double[] data){

        try{
            File file = new File(filename+".txt");
            if(!file.exists()) file.createNewFile();

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            String headerString = "#";
            String dataString = "";

            for(int i = 0; i < data.length; i++){
                bw.write(String.valueOf(data[i]));
                bw.newLine();
            }
            bw.close();

        }catch (IOException e){}

    }
}
