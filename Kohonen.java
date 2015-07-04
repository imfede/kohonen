import  java.util.*;
import java.io.*;



class Kohonen {

    static int inputDimension = 3;
    static int numberOfNeurons = 900;
    static int netDimension = 2;
    static int numberOfRows  = 30; // assert (nON % nOR == 0)
    static int numberOfSamples = 10000;
    static int maxIterations = 1000000;
    static double stopParameter = 0.0005;
    static int printStep = 20;

    static int what = 5; // 1 uniforme, 2 cerchio, 3 strana, 4 sfera vuota, 5 gaussian/sin+cos

    public static Data pickRandomSample(Data[] samples) {
	int number = (int)(Math.random() * samples.length);
	return samples[ number ];
    }

    public static Data[] createNewSampleList(int n, int dim) {
	Data[] ret = new Data[n];
	for(int i=0; i<n; i++) {
	    if( dim == 1 ) {
		double x = 0;
		
		// quadratica [-5,5]
		if( what == 3 ) {
		    double tmp = Math.random()*2 - 1;
		    x = tmp*tmp*5;
		}

		// uniforme [-5,5]
		if( what == 1 ) {
		    x = Math.random()*2 - 1;
		}

		ret[i] = new Data( x );
	    } else if( dim ==2 ) {
		double x = 0;
		double y = 0;

		// uniforme in [-5, 5]
		if( what == 1 ) {
		    //while( Math.abs(x) < 2.5 && Math.abs(y) < 2.5 ) {
			x = Math.random()*10 - 5;
			y = Math.random()*10 - 5;
			//}
		}
		
		// triangolo [-5,5]
		if( what == 3 ) {
		    x = 5;
		    y = 5;
		    while ( ((y <= -4) || ( (x+4)<=(y+4)*0.5 ) || (-x+4)<=(y+4)*0.5) ) {
			x = Math.random()*10 - 5;
			y = Math.random()*10 - 5;
		    }
		}

		// cerchio [-5,5]
		if( what == 2 ) {
		    x=5;
		    y=5;
		    while ( (x*x + y*y > 25) || (x*x + y*y < 4)) {
			x = (Math.random()-0.5)*10;
			y = (Math.random()-0.5)*10;
		    }
		}

		// cerchio vuoto [-6, 6]
		if( what == 4 ) {
		    x = Math.random()*10 - 5;
		    y = Math.signum( Math.random() - 0.5 ) * Math.sqrt( 25 - x*x );

		    Random rnd = new Random();
		    x = x + rnd.nextGaussian()*0.3;
		    y = y + rnd.nextGaussian()*0.3;
		}

		// gaussian square [-5, 5]
		if( what == 5 ) {
		    x = 5;
		    y = 5;
		    Random rnd = new Random();
		    while( Math.abs(x) >= 5 || Math.abs(y) >= 5 ) {
			x = rnd.nextGaussian();
			y = rnd.nextGaussian();			 
		    }

		}
		
		ret[i] = new Data( x, y );
	    } else if( dim == 3 ) {
		double x = 0;
		double y = 0;
		double z = 0;
		
		// uniforme [-5, 5]
		if ( what == 1 ) {
		    x = Math.random()*10 - 5;
		    y = Math.random()*10 - 5;
		    z = Math.random()*10 - 5;
		}

		// sfera piena [-5, 5]
		if( what == 2 ) {
		    x = 5;
		    y = 5;
		    z = 5;
		    while ( x*x + y*y + z*z > 25 ) {
			x = (Math.random()-0.5)*10;
			y = (Math.random()-0.5)*10;
			z = (Math.random()-0.5)*10;
		    }
		}
		
		// paraboloide [-5, 5]
		if( what == 3 ) {
		    x = Math.random()*10 - 5;
		    y = Math.random()*10 - 5;
		    z = 0.2 * (x*x + y*y) - 5;

		    // Random rnd = new Random();
		    // z = z + rnd.nextGaussian()*0.4;
		}

		// semi - sfera vuota [-5, 5]
		if( what == 4 ) {
		    x = Math.random()*10 - 5;
		    y = Math.random()*10 - 5;
		    z = -Math.sqrt( 25 - x*x - y*y );  // Math.signum( Math.random()-0.5 )*

		    //Random rnd = new Random();
		    //z = z + rnd.nextGaussian()*0.4;
		    
		    
		}

		// sin+cos [-5, 5]
		if( what == 5 ) {
		    Random rnd = new Random();
		    x = Math.random()*10 -5;
		    y = Math.random()*10 -5;
		    z = 2 * ( Math.sin(x) + Math.cos(y) ) + (rnd.nextGaussian()*0.2);
		}
		
		ret[i] = new Data( x, y, z );
	    }
	}
	return ret;
    }

    public static double fnoise(double x) {
	double noiseFactor = 0.5;
	Random rnd = new Random();
	
	// quadratica
	double a = -0.25;	
	double b = 0.75;
	double c = 1.5;
	double d = -2;
	
	double y = (a * x * x * x) + (b * x * x) + (c * x) + d;
	return y + ( rnd.nextGaussian() * noiseFactor );
    }

    public static double distance( double[] node, Data sample) {
	double sum=0;
	for(int i=0; i<node.length; i++) {
	    sum += (node[i]-sample.array[i])*(node[i]-sample.array[i]);
	}
	return Math.sqrt( sum );
    }

    public static double nodeDistance(int i, int j) {
	if( netDimension == 1 ) {
	    return (i-j);
	}
	if( netDimension == 2 ) {
	    int xi = i % numberOfRows;
	    int yi = i / numberOfRows;

	    int xj = j % numberOfRows;
	    int yj = j / numberOfRows;

	    return Math.abs(xj-xi) + Math.abs(yj-yi);
	}

	return 0;
    }
    
    public static double train( Data sample, double sigma, double[][] net, double stopParameter, double iterationCount ) {
	
	// get the nearest node
	int winnerNode = -1;
	double minDistance = Double.MAX_VALUE;
	for( int i=0; i < net.length; i++ ) {
	    double d = distance( net[i], sample );
	    if( d <= minDistance ) {
		winnerNode = i;
		minDistance = d;
	    }
	}

	double alphaCoefficent = Math.exp( - ( Math.pow(iterationCount / 2000, 2) ) ); // e^(- (x/2000)^2 )
	
	// update every node
	double maxDelta = 0;
	for( int i=0; i<net.length; i++ ) {
	    for( int j=0; j<net[0].length; j++ ) {
		double sigmaCoefficent = Math.exp( - ( Math.pow( nodeDistance(i, winnerNode), 2)  ) / ( 2*sigma*sigma) );
		double newIJ = alphaCoefficent * sigmaCoefficent * ( sample.array[j] - net[i][j] );
		if( Math.abs(newIJ) > maxDelta ) {
		    maxDelta = Math.abs(newIJ);
		}
		if( Math.abs(newIJ) > stopParameter ) {
		    net[i][j] += newIJ;
		}
	    }
	}

	return maxDelta;
    }

    public static double reduceSigma( int iterationCount, int numberOfNeurons ) {
	return numberOfNeurons * 0.20 * 1 / ( 1 + Math.pow( iterationCount * 0.01, 2 ) );
    }

    public static void printNet( double[][] net, PrintWriter writer ) {	
	if( netDimension == 1 ) {
	    for(int i=0; i<net.length; i++) {
		for(int j=0; j<net[0].length; j++) {
		    writer.print( net[i][j] + ", " );
		}
		writer.print( "0 \n" );
	    }
	} else if( netDimension == 2 ) {
	    // print rows
	    for(int i=0; i< numberOfRows; i++) {
		
		// i = number of the selected row
		for(int j = 0; j < net.length; j++) {
		    if( j % numberOfRows == i ) { // i'm on an element of the i-th row
			for(int x=0; x<net[0].length; x++) {
			    writer.print( net[j][x] + ", " );
			}
			writer.print( "0 \n");
		    }
		}
	    }

	    // print columns
	    for(int i=0; i < net.length / numberOfRows; i++) {

		// i = selected column
		for(int j = 0; j< net.length; j++) {
		    if( j / numberOfRows == i ) { //i'n on the i-th column
			for(int x=0; x<net[0].length; x++) {
			    writer.print( net[j][x] + ", " );
			}
			writer.print( "0, \n" );
		    }
		}
	    }
	}
    }

    public static void printSamples( Data[] samples, PrintWriter writer ) {
	for(int i=0; i<samples.length; i++) {
	    for( int j=0; j<samples[i].array.length; j++ ) {
	        writer.print( samples[i].array[j] + ", ");
	    }
	    writer.print("0 \n");
	}
    }
    
    public static void main(String[] args) {

	// net creation
	double[][] net = new double[ numberOfNeurons ][];
	for(int i=0; i < numberOfNeurons; i++) {
	    net[i] = new double[ inputDimension ];
	}

	// net initialisation
	if( netDimension == 1 || netDimension == 2 ) {
	    for(int i=0; i < net.length; i++) {
		for( int j=0; j < net[0].length; j++) {
		    net[i][j] = Math.random();
		}
	    }
	}
	
	
	Data[] samples = createNewSampleList( numberOfSamples, inputDimension );

	try {
	    PrintWriter writer = new PrintWriter("datas/samples", "UTF-8");
	    printSamples( samples, writer );
	    writer.close();
	}
	catch( FileNotFoundException e ) {}
	catch( UnsupportedEncodingException e) {}
	
	int iterationCount = 0;
	double sigma = reduceSigma( iterationCount, numberOfNeurons);
	double maxDelta = 1;
	int midPrint = 0;

	while( iterationCount < maxIterations && maxDelta > stopParameter ) {

	    if( iterationCount % printStep == 0 ) {
		try {
		    String file = "datas/data" + String.format( "%05d", midPrint );
		    midPrint += 1;
		    PrintWriter writer = new PrintWriter(file, "UTF-8");
		    printNet( net, writer );
		    writer.close();
		}
		catch( FileNotFoundException e ) {}
		catch( UnsupportedEncodingException e) {}
	    }
	    
	    Data sample = pickRandomSample( samples );
	    maxDelta = train( sample, sigma, net, stopParameter, iterationCount );
	    sigma = reduceSigma( iterationCount, numberOfNeurons );
	    	    
	    iterationCount += 1;
	    
	}

	try {
	    String file = "datas/data" + String.format( "%05d", midPrint );;
	    PrintWriter writer = new PrintWriter(file, "UTF-8");
	    printNet( net, writer );
	    writer.close();
	}
	catch( FileNotFoundException e ) {}
	catch( UnsupportedEncodingException e) {}
	
	System.out.println("Iterations: " + iterationCount );
    }
}

class Data {

    public double[] array;

    public Data( double x ) {
	this.array = new double[1];
	array[0]=x;
    }

    public Data( double x, double y ) {
	this.array = new double[2];
	array[0]=x;
	array[1]=y;
    }

    public Data( double x, double y, double z) {
	this.array = new double[3];
	array[0]=x;
	array[1]=y;
	array[2]=z;
    }
}
