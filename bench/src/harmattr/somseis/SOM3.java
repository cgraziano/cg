package somseis;

import java.util.Random;
import static java.lang.Math.*;

import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.SimplePlot;

public class SOM3 {
	private int iter, sIter, ic = 0, jc = 0, kc = 0;
	private float alpha = 0, radius, rStart, sigma = 0;
	private static int n_att;
	private float[][][][] nodes;
	private static float[][][] trainData;
	private float[][][][][][] distanceTable;
	float[] pickedData;
	private int n1som, n2som, n3som;

	public SOM3(int numIterations, int n1, int n2, int n3, int numAttributes) {
		this.iter = numIterations;
		this.n_att = numAttributes;

		n1som = n1;
		n2som = n2;
		n3som = n3;

		nodes = new float[n3som][n2som][n1som][numAttributes];
		randomizeAllNodesWeights(nodes);

		int max = n1;
		if (max < n2) {
			max = n2;
		}
		if (max < n3) {
			max = n3;
		}
		this.radius = max * .5f + 1.f;
		this.rStart = this.radius;

		distanceTable = new float[n3som][n2som][n1som][n3som][n2som][n1som];
		fillDistanceTable();
	}

	private static void randomizeAllNodesWeights(float[][][][] nodes) {

		int n1 = nodes.length;
		int n2 = nodes[0].length;
		int n3 = nodes[0][0].length;
		for (int i = 0; i < n1; ++i)
			for (int j = 0; j < n2; ++j)
				for (int k = 0; k < n3; ++k)
					randomizeNodeWeight(nodes[i][j][k]);
	}

	private static void randomizeNodeWeight(float[] node) {
		Random r = new Random();
		for (int i = 0; i < n_att; ++i) {
			node[i] = r.nextFloat();// .5f;
		}
	}

	private int i1(int node) {
		return node % n1som;
	}

	private int i2(int node) {
		return node / n1som;
	}

	/*
	 * Should this be private? If so how can we do that? Refer to Run.java line
	 * 132
	 */
	public int node(int i1, int i2, int i3) {
		return i1 + i2 * n1som + i3 * n2som * n1som;
	}

	private static float[] randomData() {
		Random r = new Random();
		int n1 = trainData[0].length;
		int n2 = trainData.length;
		int i1 = r.nextInt(n1);
		int i2 = r.nextInt(n2);
		return trainData[i2][i1];
	}

	public void train(float[][][] trainData) {
		this.trainData = trainData;
		iterationManager();
		System.out.println("rStart = " + rStart);
	}

	public void iterationManager() {
		for (int t = 0; t < iter; ++t) {
			/*********** Start Plotting Purposes **************/
			if (t == 0 || t == 20 || t == 100 || t == 1000 || t == 5000
					|| t == 100000) {
				float[][][] weightn11 = new float[n1som][n2som][n3som];
				float[][][] weightn12 = new float[n1som][n2som][n3som];
				float[][][] weightn21 = new float[n2som][n1som][n3som];
				float[][][] weightn22 = new float[n2som][n1som][n3som];
				float[][][] weightn31 = new float[n3som][n2som][n1som];
				float[][][] weightn32 = new float[n3som][n2som][n1som];

				for (int i = 0; i < n3som; ++i) {
					for (int j = 0; j < n2som; ++j) {
						for (int k = 0; k < n1som; ++k) {
							weightn31[i][j][k] = nodes[i][j][k][0];
							weightn32[i][j][k] = nodes[i][j][k][1];
							weightn21[j][i][k] = nodes[i][j][k][0];
							weightn22[j][i][k] = nodes[i][j][k][1];
							weightn11[k][j][i] = nodes[i][j][k][0];
							weightn12[k][j][i] = nodes[i][j][k][1];
						}
					}
				}

				String iterations = Integer.toString(t);
				PointsView pointsn1 = new PointsView(weightn11[0], weightn12[0]);
				PointsView pointsn2 = new PointsView(weightn21[0], weightn22[0]);
				PointsView pointsn3 = new PointsView(weightn31[0], weightn32[0]);

				pointsn1.setMarkColor(java.awt.Color.RED);
				pointsn1.setMarkStyle(PointsView.Mark.CROSS);
				pointsn1.setMarkSize(4);

				pointsn2.setMarkColor(java.awt.Color.RED);
				pointsn2.setMarkStyle(PointsView.Mark.CROSS);
				pointsn2.setMarkSize(4);

				pointsn3.setMarkColor(java.awt.Color.RED);
				pointsn3.setMarkStyle(PointsView.Mark.CROSS);
				pointsn3.setMarkSize(4);

				PlotPanel ppn1n2 = new PlotPanel();
				ppn1n2.addTiledView(pointsn1);
				ppn1n2.addTiledView(pointsn2);

				PlotPanel ppn2n3 = new PlotPanel();
				ppn2n3.addTiledView(pointsn2);
				ppn2n3.addTiledView(pointsn3);

				PlotPanel ppn1n3 = new PlotPanel();
				ppn1n3.addTiledView(pointsn1);
				ppn1n3.addTiledView(pointsn3);

				ppn1n2.setLimits(0, 0, 7, 7); // change to 5,5 for seismic data
				ppn2n3.setLimits(0, 0, 7, 7); // change to 5,5 for seismic data
				ppn1n3.setLimits(0, 0, 7, 7); // change to 5,5 for seismic data

				ppn1n2.addTitle(iterations);
				ppn2n3.addTitle(iterations);
				ppn1n3.addTitle(iterations);
				// pp.setLimits(0, 0, 1, 1);
				PlotFrame pfn1n2 = new PlotFrame(ppn1n2);
				PlotFrame pfn2n3 = new PlotFrame(ppn2n3);
				PlotFrame pfn1n3 = new PlotFrame(ppn1n3);

				pfn1n2.setVisible(true);
				pfn2n3.setVisible(true);
				pfn1n3.setVisible(true);
			}
			/*********** End Plotting Purposes **************/

			pickedData = randomData();
			findWinner(pickedData);
			if (t < iter * .1f) {
				alpha = (float) (0.9 * (1. - (t / (float) (iter * .1))));
				sigma = (float) (.33355f)
						* (rStart - ((rStart - 1) / (iter * .1f)) * t);// (2-(5.*t)/(float)(3.*iter));//(float)
																		// ((sRadius/3.0-1)/(-(float)iter))*((1-sRadius)/(-sRadius+3))*t+(sRadius/3f);//(2-(5.*t)/(float)(3.*iter));

			} else {
				alpha = .01f;
				sigma = .33355f;
			}

			radius = 3.f * sigma;// (float)
									// (exp(1.5*(1-t/(float)iter)));//(float)3.*sigma;//(exp(1.5*(1-t/(float)iter)));

			// System.out.println("alpha= "+alpha);
			updateNodes(ic, jc, kc);

			// System.out.println("alpha= "+alpha);

		}
	}

	public void updateNodes(int ic, int jc, int kc) {
		float h = 0;
		for (int i = 0; i < n3som; ++i) {
			for (int j = 0; j < n2som; ++j) {
				for (int k = 0; k < n1som; ++k) {
					updateNode(ic, jc, kc, i, j, k);
				}
			}
		}
	}

	public void updateNode(int ic, int jc, int kc, int ii, int jj, int kk) {
		float h = 0;
		float distance = distanceTable[ic][jc][kc][ii][jj][kk];
		if (distance == 0) {
			for (int i = 0; i < n_att; ++i) {
				// h = (float)
				// (alpha*Math.exp(-Math.pow(distance(ic,jc,ii,jj),2)/Math.pow(alpha,2)));
				h = (float) (alpha * exp(-pow(
						distanceTable[ic][jc][kc][ii][jj][kk], 2)
						/ (2. * pow(sigma, 2))));

				// System.out.println("h= "+h);
				nodes[ic][jc][kc][i] = nodes[ic][jc][kc][i] + h
						* (pickedData[i] - nodes[ic][jc][kc][i]);
			}
		} else if (distance <= radius) {
			for (int i = 0; i < n_att; ++i) {
				// h = (float)
				// (alpha*Math.exp(-Math.pow(distance(ic,jc,ii,jj),2)/Math.pow(alpha,2)));
				h = (float) (alpha * exp(-pow(
						distanceTable[ic][jc][kc][ii][jj][kk], 2)
						/ (2. * pow(sigma, 2))));
				// System.out.println("h= "+distance(ic,jc,ii,jj));
				nodes[ii][jj][kk][i] = nodes[ii][jj][kk][i] + h
						* (pickedData[i] - nodes[ii][jj][kk][i]);

			}
		}
	}

	public void train(float[][][][] trainData) {
		// Iteration Manager
		// Update
	}

	public float[][] classifyAllData(float[][][] inputData) {
		int n1 = inputData.length;
		int n2 = inputData[0].length;
		float[][] classifiedData = new float[n1][n2];

		int winner = 0;// initialize winner node
		for (int i = 0; i < n1; ++i) {
			for (int j = 0; j < n2; ++j) {
				classifiedData[i][j] = classifyOneData(inputData[i][j]);
			}
		}

		return classifiedData;
	}

	/*
	 * Changed this method to have winCatNum equal the value that is spit out
	 * from the node() method.
	 */
	private int classifyOneData(float[] inputData) {
		int n1 = nodes.length;
		int n2 = nodes[0].length;
		int n3 = nodes[0][0].length;
		int winCatNum = 0;
		float minEucDis = 999999999;
		float eucDis = 0;
		for (int i = 0; i < n1; ++i) {
			for (int j = 0; j < n2; ++j) {
				for (int k = 0; k < n3; ++k) {
					eucDis = euclideanDistance(inputData, nodes[i][j][k]);
					if (eucDis < minEucDis) {
						minEucDis = eucDis;
						winCatNum = node(k, j, i);
					}
				}
			}
		}
		System.out.println("Wincatnum= " + winCatNum);
		return winCatNum;
	}

	private void fillDistanceTable() {
		for (int i1 = 0; i1 < n3som; ++i1)
			for (int i2 = 0; i2 < n2som; ++i2)
				for (int i3 = 0; i3 < n1som; ++i3)
					for (int i4 = 0; i4 < n3som; ++i4)
						for (int i5 = 0; i5 < n2som; ++i5)
							for (int i6 = 0; i6 < n1som; ++i6) {
								distanceTable[i1][i2][i3][i4][i5][i6] = (float) (sqrt(pow(
										(i1 - i4), 2) + pow((i2 - i5), 2)) + pow(
										(i3 - i6), 2));
							}
	}

	private void findWinner(float[] singleTrainData) {
		float xmin = 99f;

		int n3 = nodes.length;
		int n2 = nodes[0].length;
		int n1 = nodes[0][0].length;
		for (int i = 0; i < n3; ++i) {
			for (int j = 0; j < n2; ++j) {
				for (int k = 0; k < n1; ++k) {
					float eucDis = euclideanDistance(singleTrainData,
							nodes[i][j][k]);
					if (eucDis < xmin) {
						xmin = eucDis;
						ic = i;
						jc = j;
						kc = k;
					}
				}
			}
		}
		// System.out.println(xmin);
	}

	private float euclideanDistance(float[] input, float[] node) {
		float eucDis;
		float diff = 0; // the euclidean distance without the square root
		for (int i = 0; i < n_att; ++i) {
			diff += (input[i] - node[i]) * (input[i] - node[i]);

		}
		eucDis = (float) Math.sqrt((double) diff);
		return eucDis;
	}

	/**
	 * The weights are re-randomized and the category numbers are re-issued.
	 */
	public void resetWeights() {
		nodes = new float[n3som][n2som][n1som][n_att];
		randomizeAllNodesWeights(nodes);
	}

	public float[][][][] getNodes() {
		return nodes;
	}
}