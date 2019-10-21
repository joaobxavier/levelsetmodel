package edu.harvard.sysbio.gavea.output;

/**
 * Operations for writing matrices to files or strings
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 26, 2007
 */
public class MatrixWriting {
	/**
	 * The data separator
	 */
	private static final String SEPARATOR = " ";

	/**
	 * Write the matrix to a string
	 * 
	 * @param matrix
	 * @return a string with the matrix (space separated values)
	 */
	public static String matrixToString(double[][] matrix) {
		StringBuffer out = new StringBuffer();
		for (int i = matrix.length - 1; i >= 0; i--) {
			for (int j = 0; j < matrix[0].length; j++) {
				out.append(matrix[i][j]);
				// change here for format (presently space separated values
				out.append(SEPARATOR);
			}
			out.append("\n");
		}
		return out.toString();
	}
}
