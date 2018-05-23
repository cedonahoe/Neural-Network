import java.util.*;
import java.io.*;
import java.io.FileNotFoundException;

public class Neural {

    public static double inA(double[] weights, double[] input) {
        double ua;
        ua = (1 * weights[0]) + (weights[1] * input[0]) + (weights[2] * input[1]);
        return ua;
    }
    public static double inB(double[] weights, double[] input) {
        double ub;
        ub = (1 * weights[3]) + (weights[5] * input[1]) + (weights[4] * input[0]);
        return ub;
    }
    public static double inC(double va, double vb, double[] weights) {
        double uc;
        uc = (1 * weights[6]) + (weights[7] * va) + (weights[8] * vb);
        return uc;
    }
    public static double outA(double in) {
        if(in > 0) {
            return in;
        }
        return 0.0;
    }
    public static double outB(double in) {
        if(in > 0) {
            return in;
        }
        return 0.0;
    }
    public static double outC(double in) {
        double c;
        double out;
        c = 1 + Math.pow(Math.E, -in);
        out = 1 / c;
        return out;

    }
    public static double Error(double v, double y) {
        double e = .5 * (Math.pow(v - y, 2));
        return e;
    }

    public static double[] Adjusted(double[] weights, double[] partials, double step) {
        double[] adj = new double[9];

        for(int i = 0; i < weights.length; i++) {
            adj[i] = weights[i] - (step * partials[i]);
        }
        return adj;
    }

    public static void main(String[] args) throws FileNotFoundException {
        int flag = Integer.valueOf(args[0]);
        double[] weights = new double[9];
        double[] adjusted = new double[9];
        double[] partials = new double[9];
        double[] inputs = new double[2];
        double yval = 0.0;
        double stepsize = 0.0;


        //print Ua Va Ub Vb with a space in between
        if(100 == flag) {
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }
            inputs[0] = Double.valueOf(args[10]);
            inputs[1] = Double.valueOf(args[11]);

            double UA = inA(weights, inputs);
            double UB = inB(weights, inputs);
            double VA = outA(UA);
            double VB = outB(UB);
            double UC = inC(VA,VB, weights);
            double VC = outC(UC);
            System.out.printf("%.5f %.5f %.5f %.5f %.5f %.5f\n", UA, VA, UB, VB, UC, VC);
        }
        //print error, partial derivative of error wrt Vc, and partial derivative of error wrt Uc
        if(200 == flag) {
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }
            inputs[0] = Double.valueOf(args[10]);
            inputs[1] = Double.valueOf(args[11]);
            yval = Double.valueOf(args[12]);

            double UA = inA(weights, inputs);
            double UB = inB(weights, inputs);
            double VA = outA(UA);
            double VB = outB(UB);
            double UC = inC(VA,VB, weights);
            double VC = outC(UC);

            double err = Error(VC, yval);
            double dvc = VC - yval;
            double duc = dvc * (VC * (1 - VC));

            System.out.printf("%.5f %.5f %.5f\n", err, dvc, duc);
        }
        //print partial derivatives of error wrt: Va, Ua, Vb, Ub separated by spaces
        if(300 == flag) {
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }
            inputs[0] = Double.valueOf(args[10]);
            inputs[1] = Double.valueOf(args[11]);
            yval = Double.valueOf(args[12]);

            double UA = inA(weights, inputs);
            double UB = inB(weights, inputs);
            double VA = outA(UA);
            double VB = outB(UB);
            double UC = inC(VA,VB, weights);
            double VC = outC(UC);

            double err = Error(VC, yval);
            double dvc = VC - yval;
            double duc = dvc * (VC * (1 - VC));

            double dva = duc * weights[7];
            double dua = 0.0;
            if(UA >= 0.0) {
                dua = 1.0 * dva;
            }
            else if(UA < 0.0) {
                dua = 0.0;
            }
            double dvb = duc * weights[8];
            double dub = 0.0;
            if(UB >= 0.0) {
                dub = 1.0 * dvb;
            }
            else if(UB < 0.0) {
                dub = 0.0;
            }
            System.out.printf("%.5f %.5f %.5f %.5f\n", dva, dua, dvb, dub);
        }
        //print the partial derivative of the error wrt each edge weight space separated
        if(400 == flag) {
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }
            inputs[0] = Double.valueOf(args[10]);
            inputs[1] = Double.valueOf(args[11]);
            yval = Double.valueOf(args[12]);

            double UA = inA(weights, inputs);
            double UB = inB(weights, inputs);
            double VA = outA(UA);
            double VB = outB(UB);
            double UC = inC(VA,VB, weights);
            double VC = outC(UC);

            double err = Error(VC, yval);
            double dvc = VC - yval;
            double duc = dvc * (VC * (1 - VC));

            double dva = duc * weights[7];
            double dua = 0.0;
            if(UA >= 0.0) {
                dua = 1.0 * dva;
            }
            else if(UA < 0.0) {
                dua = 0.0;
            }
            double dvb = duc * weights[8];
            double dub = 0.0;
            if(UB >= 0.0) {
                dub = 1.0 * dvb;
            }
            else if(UB < 0.0) {
                dub = 0.0;
            }

            partials[0] = 1.0 * dua;
            partials[1] = dua * inputs[0];
            partials[2] = dua * inputs[1];
            partials[3] = 1.0 * dub;
            partials[4] = dub * inputs[0];
            partials[5] = dub * inputs[1];
            partials[6] = 1.0 * duc;

            partials[7] = VA * duc;
            partials[8] = VB * duc;

            for(int i = 0; i < partials.length; i ++) {
                System.out.printf("%.5f ", partials[i]);
            }
            System.out.printf("\n");
        }
        //Here we perform one step of stochastic gradient descent with step size n
        //print The old weights on one line, the error under the old weights, the updated weights, and the error with the updated weights
        if(500 == flag) {
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }

            inputs[0] = Double.valueOf(args[10]);
            inputs[1] = Double.valueOf(args[11]);
            yval = Double.valueOf(args[12]);
            stepsize = Double.valueOf(args[13]);

            double UA = inA(weights, inputs);
            double UB = inB(weights, inputs);
            double VA = outA(UA);
            double VB = outB(UB);
            double UC = inC(VA,VB, weights);
            double VC = outC(UC);

            double err = Error(VC, yval);

            //print old weights
            for(int i = 0; i < weights.length; i ++) {
                System.out.printf("%.5f ", weights[i]);
            }
            System.out.printf("\n");
            //print old error
            System.out.printf("%.5f\n", err);

            //set up the array of partial derivatives wrt weights
            double dvc = VC - yval;
            double duc = dvc * (VC * (1 - VC));

            double dva = duc * weights[7];
            double dua = 0.0;
            if(UA >= 0.0) {
                dua = 1.0 * dva;
            }
            else if(UA < 0.0) {
                dua = 0.0;
            }
            double dvb = duc * weights[8];
            double dub = 0.0;
            if(UB >= 0.0) {
                dub = 1.0 * dvb;
            }
            else if(UB < 0.0) {
                dub = 0.0;
            }

            partials[0] = 1.0 * dua;
            partials[1] = dua * inputs[0];
            partials[2] = dua * inputs[1];
            partials[3] = 1.0 * dub;
            partials[4] = dub * inputs[0];
            partials[5] = dub * inputs[1];
            partials[6] = 1.0 * duc;

            partials[7] = VA * duc;
            partials[8] = VB * duc;

            //print adjusted weights
            adjusted = Adjusted(weights, partials, stepsize);
            for(int i = 0; i < adjusted.length; i++) {
                System.out.printf("%.5f ", adjusted[i]);
            }
            System.out.printf("\n");

            double newUA = inA(adjusted, inputs);
            double newUB = inB(adjusted, inputs);
            double newVA = outA(newUA);
            double newVB = outB(newUB);
            double newUC = inC(newVA,newVB, adjusted);
            double newVC = outC(newUC);

            double newE = Error( newVC,yval);

            System.out.printf("%.5f\n", newE);

        }
        //Use the neural network to predict if a student will get an A in this course given x1 is their score on hw2 and x2 is their midterm score
        //these scores are scaled to be between 0 and 1 already
        //we will run 1 epoch going over the 67 training items once
        //for each training item print x1, x2, y
        //then the updated weights
        //then the evaluation set error under the update weights
        //201 lines total will be printed
        if(600 == flag) {
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }
            stepsize = Double.valueOf(args[10]);
            File input = new File("hw2_midterm_A_train.txt");
            Scanner in = new Scanner(input);

            double[] w = weights;
            for(int j = 1; j <= 67; j++) {
                inputs[0] = in.nextDouble();
                inputs[1] = in.nextDouble();
                yval = in.nextDouble();

                double UA = inA(w, inputs);
                double UB = inB(w, inputs);
                double VA = outA(UA);
                double VB = outB(UB);
                double UC = inC(VA,VB, w);
                double VC = outC(UC);

                double Ec = Error(VC, yval);
                double EdVc = VC - yval;
                double EdUc = EdVc * (VC * (1 - VC));

                double EdVa = (EdUc * w[7]);
                double EdUa =0;
                if(UA >= 0.0) {
                    EdUa = 1.0 * EdVa;
                }
                else if(UA < 0.0) {
                    EdUa = 0.0;
                }
                double EdVb = (EdUc * w[8]);
                double EdUb = 0;
                if(UB >= 0.0) {
                    EdUb = 1.0 * EdVb;
                }
                else if (UB < 0.0) {
                    EdUb = 0.0;
                }

                w[0] = w[0] - (stepsize * (1.0 * EdUa));
                w[1] = w[1] - (stepsize * inputs[0] * EdUa);
                w[2] = w[2] - (stepsize * (inputs[1] * EdUa));
                w[3] = w[3] - (stepsize * (1.0 * EdUb));
                w[4] = w[4] - (stepsize * (inputs[0] * EdUb));
                w[5] = w[5] - (stepsize * (inputs[1] * EdUb));
                w[6] = w[6] - (stepsize * (1.0 * EdUc));
                w[7] = w[7] - (stepsize * (VA * EdUc));
                w[8] = w[8] - (stepsize * (VB * EdUc));
                System.out.printf("%.5f %.5f %.5f\n", inputs[0], inputs[1], yval);
                for(int x = 0; x < w.length; x++) {
                    System.out.printf("%.5f ", w[x]);
                }
                System.out.printf("\n");

                File eval = new File("hw2_midterm_A_eval.txt");
                Scanner ineval = new Scanner(eval);

                double eseterr = 0.0;
                double yeval;
                double[] evalin = new double[2];
                double tmp;
                for(int y = 0; y < 25; y++) {
                    evalin[0] = ineval.nextDouble();
                    evalin[1] = ineval.nextDouble();
                    yeval = ineval.nextDouble();
                    ineval.nextLine();
                    //need to calculate evaluation set error with new weights inputs and yval
                    double eUA = inA(w, evalin);
                    double eUB = inB(w, evalin);
                    double eVA = outA(eUA);
                    double eVB = outB(eUB);
                    double eUC = inC(eVA,eVB, w);
                    double eVC = outC(eUC);
                    tmp = .5 * (Math.pow((eVC - yeval), 2));
                    eseterr += tmp;
                }
                System.out.printf("%.5f\n", eseterr);
            }
        }

        //The same as 600 but we will specify a number of epochs to run the network
        if(700 == flag) {
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }
            stepsize = Double.valueOf(args[10]);
            double numepoch = Double.valueOf(args[11]);

            double[] w = weights;
            for(int k = 0; k < numepoch; k++) {

                File input = new File("hw2_midterm_A_train.txt");
                Scanner in = new Scanner(input);

                for (int j = 1; j <= 67; j++) {
                    inputs[0] = in.nextDouble();
                    inputs[1] = in.nextDouble();
                    yval = in.nextDouble();

                    double UA = inA(w, inputs);
                    double UB = inB(w, inputs);
                    double VA = outA(UA);
                    double VB = outB(UB);
                    double UC = inC(VA, VB, w);
                    double VC = outC(UC);

                    double Ec = Error(VC, yval);
                    double EdVc = VC - yval;
                    double EdUc = EdVc * (VC * (1 - VC));

                    double EdVa = (EdUc * w[7]);
                    double EdUa = 0;
                    if (UA >= 0.0) {
                        EdUa = 1.0 * EdVa;
                    } else if (UA < 0.0) {
                        EdUa = 0.0;
                    }
                    double EdVb = (EdUc * w[8]);
                    double EdUb = 0;
                    if (UB >= 0.0) {
                        EdUb = 1.0 * EdVb;
                    } else if (UB < 0.0) {
                        EdUb = 0.0;
                    }

                    w[0] = w[0] - (stepsize * (1.0 * EdUa));
                    w[1] = w[1] - (stepsize * inputs[0] * EdUa);
                    w[2] = w[2] - (stepsize * (inputs[1] * EdUa));
                    w[3] = w[3] - (stepsize * (1.0 * EdUb));
                    w[4] = w[4] - (stepsize * (inputs[0] * EdUb));
                    w[5] = w[5] - (stepsize * (inputs[1] * EdUb));
                    w[6] = w[6] - (stepsize * (1.0 * EdUc));
                    w[7] = w[7] - (stepsize * (VA * EdUc));
                    w[8] = w[8] - (stepsize * (VB * EdUc));

                }
                in.reset();
                //first print the weights at the end of the epoch
                for (int x = 0; x < w.length; x++) {
                    System.out.printf("%.5f ", w[x]);
                }
                System.out.printf("\n");

                //then print the evaluation set error once for each epoch
                File eval = new File("hw2_midterm_A_eval.txt");
                Scanner ineval = new Scanner(eval);
                double eseterr = 0.0;
                double yeval;
                double[] evalin = new double[2];
                double tmp;
                for (int y = 0; y < 25; y++) {
                    evalin[0] = ineval.nextDouble();
                    evalin[1] = ineval.nextDouble();
                    yeval = ineval.nextDouble();
                    ineval.nextLine();
                    //need to calculate evaluation set error with new weights inputs and yval
                    double eUA = inA(w, evalin);
                    double eUB = inB(w, evalin);
                    double eVA = outA(eUA);
                    double eVB = outB(eUB);
                    double eUC = inC(eVA, eVB, w);
                    double eVC = outC(eUC);
                    tmp = .5 * (Math.pow((eVC - yeval), 2));
                    eseterr += tmp;
                }
                System.out.printf("%.5f\n", eseterr);
            }
        }
        //Implemented the same as 700; we can set the number of epochs for the network to run
        //this implementation will stop running when the evaluation set error increases over one epoch
        //This will result in a well trained neural network
        //print the number of epochs that were executed. Then print the final weights, and finally the evaluation set error
        //finally print the test set classification accuracy
        if(800 == flag) {
            double lowesterr = Double.MAX_VALUE;
            for(int i = 1; i < 10; i ++) {
                weights[i - 1] = Double.valueOf(args[i]);
            }
            stepsize = Double.valueOf(args[10]);
            double numepoch = Double.valueOf(args[11]);

            double[] w = weights;
            for(int k = 0; k < numepoch; k++) {

                File input = new File("hw2_midterm_A_train.txt");
                Scanner in = new Scanner(input);

                for (int j = 1; j <= 67; j++) {
                    inputs[0] = in.nextDouble();
                    inputs[1] = in.nextDouble();
                    yval = in.nextDouble();

                    double UA = inA(w, inputs);
                    double UB = inB(w, inputs);
                    double VA = outA(UA);
                    double VB = outB(UB);
                    double UC = inC(VA, VB, w);
                    double VC = outC(UC);

                    double Ec = Error(VC, yval);
                    double EdVc = VC - yval;
                    double EdUc = EdVc * (VC * (1 - VC));

                    double EdVa = (EdUc * w[7]);
                    double EdUa = 0;
                    if (UA >= 0.0) {
                        EdUa = 1.0 * EdVa;
                    } else if (UA < 0.0) {
                        EdUa = 0.0;
                    }
                    double EdVb = (EdUc * w[8]);
                    double EdUb = 0;
                    if (UB >= 0.0) {
                        EdUb = 1.0 * EdVb;
                    } else if (UB < 0.0) {
                        EdUb = 0.0;
                    }

                    w[0] = w[0] - (stepsize * (1.0 * EdUa));
                    w[1] = w[1] - (stepsize * inputs[0] * EdUa);
                    w[2] = w[2] - (stepsize * (inputs[1] * EdUa));
                    w[3] = w[3] - (stepsize * (1.0 * EdUb));
                    w[4] = w[4] - (stepsize * (inputs[0] * EdUb));
                    w[5] = w[5] - (stepsize * (inputs[1] * EdUb));
                    w[6] = w[6] - (stepsize * (1.0 * EdUc));
                    w[7] = w[7] - (stepsize * (VA * EdUc));
                    w[8] = w[8] - (stepsize * (VB * EdUc));

                }
                in.reset();

                //need to calculate the evaluation set error and compare with the previous epochs value
                File eval = new File("hw2_midterm_A_eval.txt");
                Scanner ineval = new Scanner(eval);
                double eseterr = 0.0;
                double yeval;
                double[] evalin = new double[2];
                double tmp;
                for (int y = 0; y < 25; y++) {
                    evalin[0] = ineval.nextDouble();
                    evalin[1] = ineval.nextDouble();
                    yeval = ineval.nextDouble();
                    ineval.nextLine();
                    //need to calculate evaluation set error with new weights inputs and yval
                    double eUA = inA(w, evalin);
                    double eUB = inB(w, evalin);
                    double eVA = outA(eUA);
                    double eVB = outB(eUB);
                    double eUC = inC(eVA, eVB, w);
                    double eVC = outC(eUC);
                    tmp = .5 * (Math.pow((eVC - yeval), 2));
                    eseterr += tmp;
                }

                //need to keep track of the smallest evaluation set error so we can stop execution when it increases
                if(eseterr < lowesterr) {
                    lowesterr = eseterr;
                }
                //the eval set err has increased we should not run the network for more epochs as it is well-trained
                if(eseterr > lowesterr) {
                    //print the number of epochs executed
                    System.out.printf("%d\n", k + 1);
                    //print the final weights
                    for (int x = 0; x < w.length; x++) {
                        System.out.printf("%.5f ", w[x]);
                    }
                    System.out.printf("\n");
                    //print the evaluation set error
                    System.out.printf("%.5f\n", eseterr);

                    //print the test classification accuracy
                    File test = new File("hw2_midterm_A_test.txt");
                    Scanner intest = new Scanner(test);
                    double ytest;
                    double[] testin = new double[2];
                    double matchcount = 0.0;

                    for (int l = 0; l < 25; l++) {
                        testin[0] = intest.nextDouble();
                        testin[1] = intest.nextDouble();
                        ytest = intest.nextDouble();
                        intest.nextLine();
                        //need to calculate evaluation set error with new weights inputs and yval
                        double tUA = inA(w, testin);
                        double tUB = inB(w, testin);
                        double tVA = outA(tUA);
                        double tVB = outB(tUB);
                        double tUC = inC(tVA, tVB, w);
                        double tVC = outC(tUC);

                        //need compare network output to .5
                        if(tVC >= .5) {
                            if(ytest == 1.0) {
                                matchcount++;
                            }
                        }
                        else if(tVC < .5) {
                            if(ytest == 0.0) {
                                matchcount++;
                            }
                        }
                    }
                    double accuracy = matchcount / 25.0;
                    System.out.printf("%.5f\n", accuracy);
                    return;
                }
            }
        }
    }

}
