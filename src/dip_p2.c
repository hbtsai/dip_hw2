#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Size 723
#define SIGMA 10
#define SNP_FREQ 128

#define max(a, b) (((a)>(b))?(a):(b))
#define min(a, b) (((a)<(b))?(a):(b))
#define max3(a, b, c) ( max(a, max(b, c))) 
#define min3(a, b, c) ( min(a, min(b, c)))
#define MAXMIN(a, b, c, d, e) max3(min3(a, b, c), min3(b, c, d), min3(c, d, e))
#define MINMAX(a, b, c, d, e) min3(max3(a, b, c), max3(b, c, d), max3(c, d, e))

int write_pgm_image(char* filename, int x_dim, int y_dim, unsigned char* image)
{
	unsigned char* y = image;
	FILE* filehandle = NULL;
	filehandle = fopen(filename, "wb");
	if (filehandle) 
	{
		fprintf(filehandle, "P5\n\n%d %d 255\n", x_dim, y_dim);
		fwrite(y, 1, x_dim * y_dim, filehandle);
		fclose(filehandle);
		return 0;
	} 
	else
	  return 1;
}

int paint_histogram(int width, int height, unsigned char* image, char* filename)
{

	unsigned char *histoimg = calloc(1, width*height);
	double histo[256]={};
	int i=0, j=0;
	for(i=0; i<256; i++)
		histo[i]=0;

	for(i=0; i<width*height; i++)
		histo[image[i]]++;

	for(i=0; i<256; i++)
	{
		histo[i]/=65536;
		histo[i]*=255*40;
	}

	for(i=0; i<256; i++)
	{
		for(j=0; j<histo[i]; j++)
		{
			if(j>255)
				break;
			histoimg[i*width+j]=238;
		}
	}

	write_pgm_image(filename, width, height, histoimg);
	free(histoimg);
	return 0;
}

int add_gaussian_noise(int sigma, int width, int height, unsigned char* image)
{
	double a=0;
	int diff=0;
	int pix=0;
	int j=0, i=0;
	for(j=0; j<width*height; j++)
	{
		a=0; 
		for(i=0; i<12; i++)
			a += rand();

		diff = (int)(sigma*(double)((double)(a/RAND_MAX)-6));
		pix = (int)image[j];
		pix+=diff;

		if(pix<0)
			pix=0;
		else if(pix>255)
			pix=255;

		image[j]=pix;

	}
	return 0;
}

int add_snp_noise(int freq, int width, int height, unsigned char* image)
{
	int *noise = calloc(sizeof(int), width*height);
	int i=0;
	for(i=0; i<width*height; i++)
		noise[i]=rand()%freq;

	for(i=0; i<width*height; i++)
	{
		if(noise[i]==0)
			image[i]=0;
		else if(noise[i]==(freq-1))
			image[i]=255;
	}
	free(noise);
	return 0;
}

int comp (const void * elem1, const void * elem2) {
    unsigned char f = *((unsigned char*)elem1);
    unsigned char s = *((unsigned char*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int remove_snp_2d(int w_size, int width, int height, unsigned char* image, unsigned char* image_r)
{
	unsigned char *array = calloc(1, w_size*w_size);
	int count=0;
	int x=0, y=0, x2=0, y2=0;
	for(x=0; x<width; x++)
	{
		for(y=0; y<height; y++)
		{
			count=0;
			for(x2=x-(w_size-1)/2; x2<=x+(w_size-1)/2; x2++)
			{
				for(y2=y-(w_size-1)/2; y2<=y+(w_size-1)/2; y2++)
				{
					if(x2<0 || y2<0 || x2>width || y2>height)
						array[count] = 0;
					else
						array[count]=image[x2*width+y2];

					count++;
				}
			}

			qsort(array, w_size*w_size, sizeof(*array), comp);
			image_r[x*width+y]=array[w_size*w_size/2];
		}
	}

	free(array);
	return 0;
}

static const int kernel9[9]={ 1, 1, 1,
							  1, 1, 1,
							  1, 1, 1 };

static const int kernel10[9]={ 1, 1, 1,
							 1, 2, 1,
							 1, 1, 1 };

static const int kernel16[9]={ 1, 2, 1,
							 2, 4, 2,
							 1, 2, 1 };

int remove_gaussian10(int w_size, int width, int height, unsigned char* image, unsigned char* image_r)
{
	int idx=0;
	int value=0;
	int pix=0;
	int x=0, y=0, x2=0, y2=0;
	for(x=0; x<width; x++)
	{
		for(y=0; y<height; y++)
		{
			idx=0;
			value=0;

			for(x2=x-(w_size-1)/2; x2<=x+(w_size-1)/2; x2++)
			{
				for(y2=y-(w_size-1)/2; y2<=y+(w_size-1)/2; y2++)
				{
					pix = (int)image[x2*width+y2];
					pix *= kernel10[idx];
					value += pix;
					idx++;
				}
			}

			image_r[x*width+y] = value/10;

		}
	}

	return 0;
}
int remove_gaussian9(int w_size, int width, int height, unsigned char* image, unsigned char* image_r)
{
	int idx=0;
	int value=0;
	int pix=0;
	int x=0, y=0, x2=0, y2=0;
	for(x=0; x<width; x++)
	{
		for(y=0; y<height; y++)
		{
			idx=0;
			value=0;

			for(x2=x-(w_size-1)/2; x2<=x+(w_size-1)/2; x2++)
			{
				for(y2=y-(w_size-1)/2; y2<=y+(w_size-1)/2; y2++)
				{
					pix = (int)image[x2*width+y2];
					pix *= kernel9[idx];
					value += pix;
					idx++;
				}
			}

			image_r[x*width+y] = value/9;

		}
	}

	return 0;
}
int remove_gaussian16(int w_size, int width, int height, unsigned char* image, unsigned char* image_r)
{
	int idx=0;
	int value=0;
	int pix=0;
	int x=0, y=0, x2=0, y2=0;
	for(x=0; x<width; x++)
	{
		for(y=0; y<height; y++)
		{
			idx=0;
			value=0;

			for(x2=x-(w_size-1)/2; x2<=x+(w_size-1)/2; x2++)
			{
				for(y2=y-(w_size-1)/2; y2<=y+(w_size-1)/2; y2++)
				{
					pix = (int)image[x2*width+y2];
					pix *= kernel16[idx];
					value += pix;
					idx++;
				}
			}

			image_r[x*width+y] = value/16;

		}
	}

	return 0;
}

float psnr(int width, int height, unsigned char* _image, unsigned char* image)
{
	float mse = 0;
	int x=0, y=0;
	for(x=0; x<width; x++)
	{
		for(y=0; y<height; y++)
		{
			mse+=pow(_image[x*width+y]-image[x*width+y], 2)/(width*height);
		}
	}

	return 10*log10((255*255)/mse);
}

void split_fisheye(int width, int height, 
		unsigned char* image, unsigned char* image_r)
{
	int i=0, j=0;
	for(i=0; i<height/2; i++)
	{
		for(j=0; j<width; j++)
		{
			image_r[i*width*2+j+width]=image[i*width+j];
		}
	}
	for(i=height/2; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			image_r[(i-height/2)*2*(width)+j]=image[i*width+j];
		}
	}
}

int main(int argc, char** argv)
{
	FILE *file = NULL;
	unsigned char Imagedata[Size*Size] = {};

	char fname[1024]={};
	if(argv[1] != NULL && strlen(argv[1])>0)
		strcpy(fname, argv[1]);
	else
	{
		fprintf(stderr, "please specify filename of raw input image.\n");
		exit(-1);
	}

	if (!(file=fopen(fname,"rb")))
	{
		fprintf(stderr, "Cannot open file!\n");
		exit(1);
	}
	fread(Imagedata, sizeof(unsigned char), Size*Size, file);
	fclose(file);

	/* save the original image for comparision */
	write_pgm_image("sample2.pgm", Size, Size, Imagedata);

	unsigned char pan1[Size*Size]={};
	split_fisheye(Size, Size, Imagedata, pan1);
	write_pgm_image("panaroma1.pgm", (2*Size), (Size/2), pan1);

	/*
	unsigned char rot[Size*Size]={};
	rotate(Size, Size, Imagedata, rot);
	write_pgm_image("rotate.pgm", Size, Size, rot);
	*/

	exit(0);
	return 0;
}


