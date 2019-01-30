import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class Toolbox {

    public static double[] averagedResults(int[][] inputData){

        int nReps = inputData.length;
        int nCounters = inputData[0].length;

        double[] averagedResults = new double[nCounters];

        //iterate over the counters, checking all the reps, then moving to next counter
        for(int c = 0; c < nCounters; c++){
            double runningTotal = 0.;
            for(int r = 0; r < nReps; r++){
                runningTotal += inputData[r][c];
            }
            averagedResults[c] = runningTotal/nReps;
        }
        return averagedResults;
    }


    public static void writeSimpleArrayAndHeadersToFile(String filename, String[] headers, double[] data){

        try{
            File file = new File(filename+".txt");
            if(!file.exists()) file.createNewFile();

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            String headerString = "#";
            String dataString = "";

            for(int i = 0; i < headers.length; i++){
                headerString+=headers[i]+"\t";
                dataString+=data[i]+"\t";
            }

            bw.write(headerString.trim());
            bw.newLine();
            bw.write(dataString.trim());
            bw.close();


        }catch (IOException e){}
    }
}
