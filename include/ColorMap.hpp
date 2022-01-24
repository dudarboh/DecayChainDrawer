#ifndef COLORMAP_H_
#define COLORMAP_H_

/**
 * This class has been obtained from 
 * https://github.com/iLCSoft/CEDViewer/blob/master/include/ColorMap.h
 * under the GNU licence.
 * Minor modifications by:
 * @author: S.Daraszewicz
 * @date: 28.08.08
 */

typedef void (*colorMapFunc)(unsigned int*,float,float,float);

typedef struct RgbColor
{
    double r;
    double g;
    double b;
} RgbColor;

typedef struct HsvColor
{
    double h;
    double s;
    double v;
} HsvColor;

class ColorMap{
 public:
  static void colorMap(unsigned int *rgb,float value,float min,float max);
  static void hotColorMap(unsigned int *rgb,float value,float min,float max);
  static void coldColorMap(unsigned int *rgb,float value,float min,float max);
  static void jetColorMap(unsigned int *rgb,float value,float min,float max);
  static void cyclicColorMap(unsigned int *rgb,float value,float min,float max);
  static void randColorMap(unsigned int *rgb,float value,float min,float max);
  static void grayColorMap(unsigned int *rgb,float value,float min,float max);
  static void blueColorMap(unsigned int *rgb,float value,float min,float max);
  static colorMapFunc selectColorMap(int cmp);
  static int RGB2HEX(int red, int green, int blue);
  static RgbColor HsvToRgb(HsvColor in);
  static unsigned long NumberToTemperature(double value, double min, double max, double s, double v);
};

#endif
