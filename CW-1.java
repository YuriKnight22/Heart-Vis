import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.lang.Math;

class Volume {
    double[][][] data;

    void Volume(int sizeX, int sizeY, int sizeZ, double scaleFactor, int chunkSize) {
        data = new double[sizeZ][sizeY][sizeX];
        for (int zChunk = 0; zChunk < sizeZ; zChunk += chunkSize) {
            int zEnd = Math.min(zChunk + chunkSize, sizeZ);
            for (int z = zChunk; z < zEnd; z++) {
                for (int y = 0; y < sizeY; y++) {
                    for (int x = 0; x < sizeX; x++) {
                        double eqResult = -Math.pow(Math.pow(x / scaleFactor, 2) + 2 * Math.pow(y / scaleFactor, 2) + Math.pow(z / scaleFactor, 2) - 1, 3)
                                        - Math.pow(x / scaleFactor, 2) * Math.pow(z / scaleFactor, 3)
                                        - 0.1 * Math.pow(y / scaleFactor, 2) * Math.pow(z / scaleFactor, 3);
                        data[z][y][x] = eqResult <= 0 ? 0 : 255;
                    }
                }
            }
        }
    }

    double[][][][] Gradient() {
        int dimX = data[0][0].length;
        int dimY = data[0].length;
        int dimZ = data.length;
        double[][][][] gradient = new double[dimZ][dimY][dimX][3];
        for (int k = 1; k < dimZ - 1; k++) {
            for (int j = 1; j < dimY - 1; j++) {
                for (int i = 1; i < dimX - 1; i++) {
                    gradient[k][j][i][0] = (data[k][j][i + 1] - data[k][j][i - 1]) / 2.0;
                    gradient[k][j][i][1] = (data[k][j + 1][i] - data[k][j - 1][i]) / 2.0;
                    gradient[k][j][i][2] = (data[k + 1][j][i] - data[k - 1][j][i]) / 2.0;
                }
            }
        }
        return gradient;
    }

    public int[][] Render(int isovalue, boolean positiveDirection, int resolution) {
        int dimX = data[0][0].length;
        int dimY = data[0].length;
        int dimZ = data.length;
        int[][] image = new int[resolution][resolution];

        for (int y = 0; y < resolution; y++) {
            for (int x = 0; x < resolution; x++) {
                // Map pixel coordinates to volume space
                double xVolume = (double)x / resolution * dimX;
                double yVolume = (double)y / resolution * dimY;

                // Ray direction along the z-axis
                double zVolume;
                if (positiveDirection)
                    zVolume = dimZ - 1;
                else
                    zVolume = 1;

                // Initialize color to black
                int color = 250;

                // Ray casting
                while (zVolume >= 0 && zVolume < dimZ) {
                    int xIndex = (int)Math.floor(xVolume);
                    int yIndex = (int)Math.floor(yVolume);
                    int zIndex = (int)Math.floor(zVolume);

                    // Check if within volume bounds
                    if (xIndex >= 0 && xIndex < dimX && yIndex >= 0 && yIndex < dimY && zIndex >= 0 && zIndex < dimZ) {
                        // Check if voxel value exceeds the isovalue
                        if (data[zIndex][yIndex][xIndex] > isovalue) {
                            // Set color based on voxel value
                            color = (int)data[zIndex][yIndex][xIndex];
                            break; // Exit loop if a voxel exceeds the isovalue
                        }
                    }

                    // Move along the ray direction
                    if (positiveDirection)
                        zVolume--;
                    else
                        zVolume++;
                    xVolume += 0.1; // Adjust the step size for better visualization
                    yVolume += 0.1;
                }

                // Assign color to the image
                image[y][x] = color;
            }
        }

        return image;
    }
}

public class CW {
    public static void saveImage(String name, int resolution, Volume volume) {
        int[][] image = new int[resolution][resolution];
        double scaleFactor = resolution * 0.30;

        for (int y = 0; y < resolution; y++) {
            for (int x = 0; x < resolution; x++) {
                double xScaled = x - resolution / 2.0;
                double yScaled = y - resolution / 2.0;

                double eqResult = -Math.pow(Math.pow(xScaled / scaleFactor, 2) + 2 * Math.pow(yScaled / scaleFactor, 2) - 1, 3)
                                    - Math.pow(xScaled / scaleFactor, 2) * Math.pow(yScaled / scaleFactor, 3)
                                    - 0.1 * Math.pow(yScaled / scaleFactor, 2) * Math.pow(Math.abs(xScaled) / scaleFactor, 3);

                if (eqResult <= 0) {
                    image[y][x] = 0;
                } else {
                    double interpolatedValue = trilinearInterpolation(volume, xScaled / scaleFactor, yScaled / scaleFactor, 0);
                    image[y][x] = (int)(interpolatedValue * 255);
                }
            }
        }
       
        BufferedImage bufferedImage = new BufferedImage(resolution, resolution, BufferedImage.TYPE_BYTE_GRAY);
        for (int y = 0; y < resolution; y++) {
            for (int x = 0; x < resolution; x++) {
                bufferedImage.setRGB(x, y, image[y][x] * 256 * 256 + image[y][x] * 256 + image[y][x]);
            }
        }

        File f = new File(name);
        try {
            ImageIO.write(bufferedImage, "tiff", f);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double trilinearInterpolation(Volume volume, double x, double y, double z) {
        int dimX = volume.data[0][0].length;
        int dimY = volume.data[0].length;
        int dimZ = volume.data.length;

        int x0 = Math.max(0, Math.min((int) Math.floor(x), dimX - 1));
        int x1 = Math.max(0, Math.min(x0 + 1, dimX - 1));
        int y0 = Math.max(0, Math.min((int) Math.floor(y), dimY - 1));
        int y1 = Math.max(0, Math.min(y0 + 1, dimY - 1));
        int z0 = Math.max(0, Math.min((int) Math.floor(z), dimZ - 1));
        int z1 = Math.max(0, Math.min(z0 + 1, dimZ - 1));

        double xd = x - x0;
        double yd = y - y0;
        double zd = z - z0;

        double c000 = volume.data[z0][y0][x0];
        double c001 = volume.data[z0][y0][x1];
        double c010 = volume.data[z0][y1][x0];
        double c011 = volume.data[z0][y1][x1];
        double c100 = volume.data[z1][y0][x0];
        double c101 = volume.data[z1][y0][x1];
        double c110 = volume.data[z1][y1][x0];
        double c111 = volume.data[z1][y1][x1];

        double c00 = c000 * (1 - xd) + c001 * xd;
        double c01 = c010 * (1 - xd) + c011 * xd;
        double c10 = c100 * (1 - xd) + c101 * xd;
        double c11 = c110 * (1 - xd) + c111 * xd;

        double c0 = c00 * (1 - yd) + c01 * yd;
        double c1 = c10 * (1 - yd) + c11 * yd;

        return c0 * (1 - zd) + c1 * zd;
    }

    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Usage: java CW <grid_resolution> <image_resolution>");
            return;
        }

        int gridResolution = Integer.parseInt(args[0]);
        int imageResolution = Integer.parseInt(args[1]);
        int chunkSize = 10; // Adjust the chunk size as needed

        Volume volume = new Volume();
        double scaleFactor = imageResolution * 0.3; // Adjust the scaling factor as needed
        volume.Volume(gridResolution, gridResolution, gridResolution, scaleFactor, chunkSize);

        // java CW 256 512
        saveImage("result.tiff", imageResolution, volume);
    }
}
