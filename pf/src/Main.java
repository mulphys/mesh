import java.io.*;

public class Main {//Read Parameters
	public static void main(String argv[]) throws IOException {
		Mesh2D mesh;
		if (argv.length != 0) {
			System.err.println("Usage: cmd < input-stream > output-stream");
			System.exit(1);
		}
		mesh = new Mesh2D();
		mesh.ReadBoundary(System.in);
		mesh.Init();
		mesh.TriangulatePF();
		mesh.Output(System.out);
		System.exit(0);
	}
}
