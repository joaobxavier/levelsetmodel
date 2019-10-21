package edu.harvard.sysbio.gavea.output;

import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;

import edu.harvard.sysbio.gavea.utils.MatrixMath;

/**
 * Keeps the path to the output directory and manages the cration of files
 * 
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 27, 2007
 */
public class OutputDirectory {
	private static final int MAX_THREADS = 7;

	private static int threadcount = 0;

	private static final String RESDIR = "results";

	private File _resultDirectory;

	/**
	 * Implements a Thread that writes the file to disk while the model
	 * continues running
	 * 
	 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 30, 2007
	 */
	private class WritingThread extends Thread {

		private double[][] _m;

		private String _name;

		boolean _copyDone = false;

		/**
		 * Copies the matrix to a new instance.
		 * 
		 * @param m
		 * @param name
		 */
		public WritingThread(double[][] m, String name) {
			int n = m.length;
			this._m = new double[n][n];
			MatrixMath.copyTo(m, _m);
			this._name = name;
			_copyDone = true;
			threadcount++;
		}

		public boolean doneCopying() {
			return _copyDone;
		}

		public void run() {
			writeMatrix(_m, _name);
			// TODO delete
			// System.out.println("finished writing " + _name);
			threadcount--;
		}
	}

	/**
	 * Create a new output directory inside directory defined by "directoryPath"
	 * 
	 * @param directoryPath
	 */
	public OutputDirectory(String directoryPath) {
		File dir = new File(directoryPath);
		try {
			if (dir.isDirectory()) {
				// append a number that increments the count of files already in
				// the subdirectory
				int n = countExistingResultsSubdirectories(dir);
				_resultDirectory = new File(directoryPath + "/" + RESDIR
						+ (n + 1));
				if (!_resultDirectory.mkdir())
					throw new RuntimeException("Error trying to create "
							+ _resultDirectory.getName());
			} else
				throw new RuntimeException("Directory " + directoryPath
						+ " not found");
			dir.mkdir();
		} catch (SecurityException e) {
			throw new RuntimeException("Error trying to write to "
					+ directoryPath);
		}
	}

	/**
	 * @param directoryPath
	 * @return the number of files or directories already existing in the path
	 */
	private int countExistingResultsSubdirectories(File directoryPath) {
		// This filter only returns text files
		FileFilter fileFilter = new FileFilter() {
			public boolean accept(File file) {
				return file.getName().contains(RESDIR);
			}
		};
		// list files in directory
		File[] files = directoryPath.listFiles(fileFilter);
		return files.length;
	}

	public void writeMatrix(double[][] m, String name) {
		File f = new java.io.File(_resultDirectory.getAbsolutePath() + "/"
				+ name);
		try {
			FileWriter fr = new FileWriter(f);
			fr.write(MatrixWriting.matrixToString(m));
			fr.close();
		} catch (IOException e) {
			throw new RuntimeException("Error writing file " + name + " in "
					+ _resultDirectory.getAbsolutePath());
		}
	}

	/**
	 * Copies the matrix to write and then creates a new thread that writes to
	 * disk while the program continues running
	 * 
	 * @param m
	 * @param name
	 */
	public void createNewThreadToWriteMatrix(double[][] m, String name) {
		//
		while (threadcount > MAX_THREADS) {
			try {
				System.out.println("waiting...");
				Thread.sleep(100);
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
		}
		//
		WritingThread writingThread = new WritingThread(m, name);
		writingThread.start();
	}
}
