#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Size 512

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
	{
	  return 1;
	}
}


static const int kernel25[25]={ 2, 4, 5, 4, 2,
							 4, 9, 12, 9, 4,
							 5, 12, 15, 12, 5,
							 4, 9, 12, 9, 4, 
							 2, 4, 5, 4, 2 };

int gaussian_filter(int w_size, int width, int height, unsigned char* image, unsigned char* image_r)
{
	int idx=0;
	int value=0;
	int pix=0;
	int x=0, y=0, x2=0, y2=0;
	for(x=0; x<height; x++)
	{
		for(y=0; y<width; y++)
		{
			idx=0;
			value=0;

			for(x2=x-(w_size-1)/2; x2<=x+(w_size-1)/2; x2++)
			{
				for(y2=y-(w_size-1)/2; y2<=y+(w_size-1)/2; y2++)
				{
					pix = (int)image[x2*width+y2];
					pix *= kernel25[idx];
					value += pix;
					idx++;
				}
			}

			image_r[x*width+y] = value/159;

		}
	}

	return 0;
}

void row_gradient(int width, int height, unsigned char* image, unsigned char* image_r, int* row_v)
{
	int i =0, j=0;
	int K=2;
	int A0=0;
	int A1=0;
	int A2=0;
	int A3=0;
	int A4=0;
	int A5=0;
	int A6=0;
	int A7=0;
	int min=0;
	for(i = 0; i<height; i++)
	{
		for(j =0; j<width; j++)
		{

			if(i==0 || j==0 || i==(width-1) || j==(height-1))
			{
				A0 = image[ i*width + j ];
				A1 = image[ i*width + j ];
				A2 = image[ i*width + j ];
				A3 = image[ i*width + j ];
				A4 = image[ i*width + j ];
				A5 = image[ i*width + j ];
				A6 = image[ i*width + j ];
				A7 = image[ i*width + j ];
			}
			else
			{
				A0 = image[ (i-1)*width + (j-1)];
				A1 = image[ (i-1)*width + j ];
				A2 = image[ (i-1)*width + (j+1) ];
				A3 = image[ (i)*width + (j+1)];
				A4 = image[ (i+1)*width + (j+1) ];
				A5 = image[ (i+1)*width + (j) ];
				A6 = image[ (i+1)*width + (j-1) ];
				A7 = image[ (i)*width + (j-1) ];
			}

			/*
			if(( ( (A0+K*A1+A2) - (A6+K*A5+A4) ) / (K+2)) < min)
				min = ( ( (A0+K*A1+A2) - (A6+K*A5+A4) ) / (K+2)) ;
			*/

			row_v[i*width+j] = (( (A2+K*A3+A4) - (A0+K*A7+A6) ) / (K+2));
			image_r[i*width+j] =(( (A2+K*A3+A4) - (A0+K*A7+A6) ) / (K+2)) + 128;
		}
	}

}

void column_gradient(int width, int height, unsigned char* image, unsigned char* image_r, int* col_v)
{
	int i =0, j=0;
	int K=2;
	int A0=0;
	int A1=0;
	int A2=0;
	int A3=0;
	int A4=0;
	int A5=0;
	int A6=0;
	int A7=0;
	int min = 0;
	for(i = 0; i<height; i++)
	{
		for(j =0; j<width; j++)
		{

			if(i==0 || j==0 || i==(width-1) || j==(height-1))
			{
				A0 = image[ i*width + j ];
				A1 = image[ i*width + j ];
				A2 = image[ i*width + j ];
				A3 = image[ i*width + j ];
				A4 = image[ i*width + j ];
				A5 = image[ i*width + j ];
				A6 = image[ i*width + j ];
				A7 = image[ i*width + j ];
			}
			else
			{
				A0 = image[ (i-1)*width + (j-1)];
				A1 = image[ (i-1)*width + j ];
				A2 = image[ (i-1)*width + (j+1) ];
				A3 = image[ (i)*width + (j+1)];
				A4 = image[ (i+1)*width + (j+1) ];
				A5 = image[ (i+1)*width + (j) ];
				A6 = image[ (i+1)*width + (j-1) ];
				A7 = image[ (i)*width + (j-1) ];
			}

			/*
			if(( ( (A0+K*A1+A2) - (A6+K*A5+A4) ) / (K+2)) < min)
				min = ( ( (A0+K*A1+A2) - (A6+K*A5+A4) ) / (K+2)) ;
			*/

			col_v[i*width+j] = ( ( (A0+K*A1+A2) - (A6+K*A5+A4) ) / (K+2)) ;
			image_r[i*width+j] = ( ( (A0+K*A1+A2) - (A6+K*A5+A4) ) / (K+2)) + 128;
		}
	}
}


void gradient(int width, int height, int* row_g, int* col_g, unsigned char* image_r)
{
	int i=0, j=0;
	for(i = 0; i<height; i++)
		for(j =0; j<width; j++)
			image_r[i*width+j] = sqrt(pow(row_g[i*width+j], 2)+ pow(col_g[i*width+j], 2));
}

void orientation(int width, int height, int* row_g, int* col_g, int* atan_g)
{
	int i=0, j=0;
	int min = 0;
	int max = 255;
	double *atan_r = calloc(sizeof(double), width*height);
	for(i = 0; i<height; i++)
		for(j =0; j<width; j++)
			atan_g[i*width+j] = (int)(atan2(col_g[i*width+j], row_g[i*width+j])/M_PI*180);

}

void supress(int width, int height, int* atan_g, unsigned char* image, unsigned char* image_r)
{
	int i=0, j=0;
	for(i =0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{


			//if(abs(atan_g[i*width+j])<6 || (180-abs(atan_g[i*width+j]))<6 )
			if(abs(atan_g[i*width+j])<22.5 || (180-abs(atan_g[i*width+j]))<22.5 )
			{
				/* compare horizontally*/
				if(j==0 || j==255)
					continue;

				if(image[i*width+j] > image[i*width+(j-1)] && image[i*width+j] > image[i*width+(j+1)] )
					image_r[i*width+j] = image[i*width+j];
				else
					image_r[i*width+j] = 0;

				continue;
			}

			//if(abs(abs(atan_g[i*width+j])-90) < 29)
			if(abs(abs(atan_g[i*width+j])-90) < 22.5)
			{
				/* compare vertically*/
				if(i==0 || i==255)
					continue;

				if(image[i*width+j] > image[(i-1)*width+j] && image[i*width+j] > image[(i+1)*width+j] )
					image_r[i*width+j] = image[i*width+j];
				else
					image_r[i*width+j] = 0;

				continue;
			}

			if(atan_g[i*width+j]>0)
			{
				//if(abs(atan_g[i*width+j]-32)<29)
				if(abs(atan_g[i*width+j]-45)<22.5)
				{
					/*compare right-up to left-down*/
					if((i==0 && j==255) || (i==255 && j==0))
						continue;

					if(image[i*width+j] > image[(i-1)*width+(j+1)] && image[i*width+j] > image[(i+1)*width+(j-1)] )
						image_r[i*width+j] = image[i*width+j];
					else
						image_r[i*width+j] = 0;

					continue;
				}

				//if(abs(atan_g[i*width+j]-148)<29)
				if(abs(atan_g[i*width+j]-135)<22.5)
				{
					/*compare left-up to right-down*/
					if((i==0 && j==0) || (i==255 && j==255))
						continue;

					if(image[i*width+j] > image[(i-1)*width+(j-1)] && image[i*width+j] > image[(i+1)*width+(j+1)] )
						image_r[i*width+j] = image[i*width+j];
					else
						image_r[i*width+j] = 0;

					continue;
				}
			}
			else
			{
				//if(abs(abs(atan_g[i*width+j])-32)<29)
				if(abs(abs(atan_g[i*width+j])-45)<22.5)
				{
					/*compare left-up to right-down*/
					if((i==0 && j==0) || (i==255 && j==255))
						continue;

					if(image[i*width+j] > image[(i-1)*width+(j-1)] && image[i*width+j] > image[(i+1)*width+(j+1)] )
						image_r[i*width+j] = image[i*width+j];
					else
						image_r[i*width+j] = 0;

					continue;
				}

			//	if(abs(abs(atan_g[i*width+j])-148)<29)
				if(abs(abs(atan_g[i*width+j])-135)<22.5)
				{
					/*compare right-up to left-down*/
					if((i==0 && j==255) || (i==255 && j==0))
						continue;

					if(image[i*width+j] > image[(i-1)*width+(j+1)] && image[i*width+j] > image[(i+1)*width+(j-1)] )
						image_r[i*width+j] = image[i*width+j];
					else
						image_r[i*width+j] = 0;

					continue;
				}
			}

			

		
		}
	}
}

void double_thresholding(int width, int height, unsigned char* image)
{
	int i=0, j=0;
	int TH=17, TL=5;
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			if(image[i*width+j]>TH)
				image[i*width+j]=255;
			else if(image[i*width+j]>TL)
				image[i*width+j]=128;
			else
				image[i*width+j]=0;
		}
	}
}

int decrease_brightness(int x, int y, unsigned char* image)
{
	int i=0; 
	for(i=0; i<x*y; i++)
		image[i]/=2;
	return 0;
}

int histogram_equalizer(int x, int y, unsigned char* image)
{
	/* create histogram */
	double histogram[256]={};
	int i=0;
	for(i=0; i<x*y; i++)
		histogram[image[i]]++;

	/* convert histogram to percentage (0-1) */
	for(i=0; i<256; i++)
		histogram[i]/=65536;

	/* cumulative histogram */
	for(i=0; i<256; i++)
		histogram[i]+=histogram[i-1];

	/* convert back to dynamic range 0-255 */
	for(i=0; i<256; i++)
		histogram[i]=(histogram[i]*254)+0.5;

	/* assign back using new histogram */
	for(i=0; i<x*y; i++)
		image[i]=histogram[image[i]];

	return 0;
}

int local_histogram_equalizer(int w_size, int width, int height, 
		int x, int y, unsigned char* image, unsigned char* image_r)
{
	// assume grayscale
	double histo[256]={};
	int i=0, x2=0, y2=0;
	for( i=0; i<256; i++)
		histo[i]=0;

	int pixels=0;
	for(x2= x-(w_size-1)/2; x2 <= x+(w_size-1)/2; x2++)
	{
		for(y2= y-(w_size-1)/2; y2 <= y+(w_size-1)/2; y2++)
		{
			if(x2<0 || y2<0 || x2>=width || y2>=height)
				continue;

			histo[image[x2*width+y2]]++;
			pixels++;
		}
	}

	/* convert histogram to percentage (0-1) */
	for(i=0; i<256; i++)
		histo[i]/=pixels;

	/* cumulative histogram */
	for(i=0; i<256; i++)
		histo[i]+=histo[i-1];

	/* convert back to dynamic range 0-255 */
	for(i=0; i<256; i++)
		histo[i]=histo[i]*255;

	/* assign back using new histogram */
	image_r[x*width+y]=histo[image[x*width+y]];
	return 0;

}

int paint_histogram(int width, int height, unsigned char* image, char* filename)
{

	unsigned char *histoimg = calloc(1, 256*256);
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
			histoimg[i*256+j]=238;
		}

		/*
		if(histo[i]>200)
		{
			fprintf(stderr, " pixel value %d is major\n", i);
		}
		*/

	}

	write_pgm_image(filename, 256, 256, histoimg);
	free(histoimg);
	histoimg=NULL;
	return 0;
}


int convert_to_black_n_white(int threshold, int width, int height, unsigned char* image)
{
	int i=0;
	for(i=0; i<width*height; i++)
	{
		if(image[i]>=threshold)
			image[i]=255;
		else
			image[i]=0;
	}
	return 0;
}

void connected(int width, int height, unsigned char* image)
{
	int i=0, j=0;
	int p0=0, p1=0 ,p2=0,p3=0,p4=0,p5=0,p6=0,p7=0;

	int changed = 1;
	while(changed)
	{
		changed = 0;
		for(i=0; i<height; i++)
		{
			for(j=0; j<width; j++)
			{
				if(i==0 || i==255 || j==0 || j==255)
					continue;
	
				p0 = image[(i-1)*width+(j-1)];
				p1 = image[(i-1)*width+(j)];
				p2 = image[(i-1)*width+(j+1)];
				p3 = image[(i)*width+(j+1)];
				p4 = image[(i+1)*width+(j+1)];
				p5 = image[(i+1)*width+(j)];
				p6 = image[(i+1)*width+(j-1)];
				p7 = image[(i)*width+(j-1)];
	
				if(image[i*width+j]==255)
				{
					if(p0==128)
					{
						changed = 1;
						image[(i-1)*width+(j-1)]=255;
					}
					if(p1==128)
					{
						changed=1;
						image[(i-1)*width+(j)]=255;
					}
					if(p2==128)
					{
						changed=1;
						image[(i-1)*width+(j+1)]=255;
					}
					if(p3==128)
					{
						changed=1;
						image[(i)*width+(j+1)]=255;
					}
					if(p4==128)
					{
						changed =1;
						image[(i+1)*width+(j+1)]=255;
					}
					if(p5==128)
					{
						changed =1;
						image[(i+1)*width+(j)]=255;
					}
					if(p6==128)
					{
						changed=1;
						image[(i+1)*width+(j-1)]=255;
					}
					if(p7==128)
					{
						changed=1;
						image[(i)*width+(j-1)]=255;
					}
					
				}
			}
		}
		if(changed)
		{
			for(i=width-1; i>0; i--)
			{
				for(j=height-1; j>0; j--)
				{
					p0 = image[(i-1)*width+(j-1)];
					p1 = image[(i-1)*width+(j)];
					p2 = image[(i-1)*width+(j+1)];
					p3 = image[(i)*width+(j+1)];
					p4 = image[(i+1)*width+(j+1)];
					p5 = image[(i+1)*width+(j)];
					p6 = image[(i+1)*width+(j-1)];
					p7 = image[(i)*width+(j-1)];

				if(image[i*width+j]==255)
				{
					if(p0==128)
					{
						changed = 1;
						image[(i-1)*width+(j-1)]=255;
					}
					if(p1==128)
					{
						changed=1;
						image[(i-1)*width+(j)]=255;
					}
					if(p2==128)
					{
						changed=1;
						image[(i-1)*width+(j+1)]=255;
					}
					if(p3==128)
					{
						changed=1;
						image[(i)*width+(j+1)]=255;
					}
					if(p4==128)
					{
						changed =1;
						image[(i+1)*width+(j+1)]=255;
					}
					if(p5==128)
					{
						changed =1;
						image[(i+1)*width+(j)]=255;
					}
					if(p6==128)
					{
						changed=1;
						image[(i+1)*width+(j-1)]=255;
					}
					if(p7==128)
					{
						changed=1;
						image[(i)*width+(j-1)]=255;
					}

				}

			}
	
			}
		}
	}

	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			if(image[i*width+j]==128)
				image[i*width+j]=0;
		}
	}
}


static const int LoG[25]={ 0, 0, -1, 0, 0,
						   0, -1, -2, -1, 0,
						   -1, -2, 16, -2, -1,
						   0, -1, -2, -1, 0,
						   0, 0, -1, 0, 0
							};

void laplacian(int w_size, int width, int height, 
		unsigned char* image, unsigned char* image_r, int* val_t)
{
	int idx=0;
	int value=0;
	int pix=0;
	int min=0;
	int max=255;
//	int val_t[Size*Size]={};
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
					pix *= LoG[idx];
					value += pix;
					idx++;
				}
			}

			if(value<min)
				min=value;

			if(value>max)
				max=value;

			val_t[x*width+y]=value;
		}
	}

	max += abs(min);

	for(x=0; x<width; x++)
		for(y=0; y<width; y++)
			image_r[x*width+y]=((val_t[x*width+y]+abs(min)))*255/max;
}

void thresh_lap(int width, int height, unsigned char* image, unsigned char* image_r, int thresh_up, int thresh_low)
{
	/* majority value of lap pixels are 175-207, mid: 186, thresh: 11*/
	int i=0, j=0;
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			if(image[i*width+j]>thresh_up)
				image_r[i*width+j]=0;
			else if(image[i*width+j]<thresh_low)
				image_r[i*width+j]=0;
			else
				image_r[i*width+j]=255;
		}
	}
}

void zero_crossing(int width, int height, int* val_t, unsigned char* image, unsigned char * image_r)
{
	int i=0, j=0;
	int p0=0;
	int p1=0;
	int p2=0;
	int p3=0;
	int p4=0;
	int p5=0;
	int p6=0;
	int p7=0;
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			if(image[i*width+j]==0)
			{
				p0 = val_t[(i-1)*width+(j-1)];
				p1 = val_t[(i-1)*width+(j)];
				p2 = val_t[(i-1)*width+(j+1)];
				p3 = val_t[(i)*width+(j+1)];
				p4 = val_t[(i+1)*width+(j+1)];
				p5 = val_t[(i+1)*width+(j)];
				p6 = val_t[(i+1)*width+(j-1)];
				p7 = val_t[(i)*width+(j-1)];

				if( ((p0<0&&p4>0) || (p0>0&&p4<0)) ||
					((p7<0&&p3>0) || (p7>0&&p3<0)) ||
					((p6<0&&p2>0) || (p6>0&&p2<0)) ||
					((p1<0&&p5>0) || (p1>0&&p5<0))
				  )
					image_r[i*width+j]=255;


			}

		}
	}
}


int main(int argc, char** argv)
{
	FILE *file = NULL;
	// image data array
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
	write_pgm_image("sample1.pgm", Size, Size, Imagedata);


	unsigned char canny_1[Size*Size] = {};
	gaussian_filter(5, Size, Size, Imagedata, canny_1);
	write_pgm_image("canny_1.pgm", Size, Size, canny_1);

	/* canny_1 is gaussian of original image */

	unsigned char canny_21[Size*Size] = {};
	int row_v[Size*Size] = {};
	int col_v[Size*Size] = {};
	row_gradient(Size, Size, canny_1, canny_21, row_v);
	write_pgm_image("canny_21.pgm", Size, Size, canny_21);

	unsigned char canny_22[Size*Size] = {};
	column_gradient(Size, Size, canny_1, canny_22, col_v);
	write_pgm_image("canny_22.pgm", Size, Size, canny_22);
	
	unsigned char canny_2[Size*Size] = {};
	gradient(Size, Size, row_v, col_v, canny_2);
	write_pgm_image("canny_2.pgm", Size, Size, canny_2);

	/*
	 * canny_2 is sobel edge detection.
	 */

	int atan[Size*Size] = {};
	orientation(Size, Size, row_v, col_v, atan);

	unsigned char canny_3[Size*Size]={};
	supress(Size, Size, atan, canny_2, canny_3);
	write_pgm_image("canny_3.pgm", Size, Size, canny_3);

	unsigned char canny_4[Size*Size]={};
	memcpy(canny_4, canny_3, sizeof(canny_3));
	double_thresholding(Size, Size, canny_4);
	write_pgm_image("canny_4.pgm", Size, Size, canny_4);

	unsigned char canny_5[Size*Size]={};
	memcpy(canny_5, canny_4, sizeof(canny_4));
	connected(Size, Size, canny_5);
	write_pgm_image("canny_5.pgm", Size, Size, canny_5);

	unsigned char lap[Size*Size]={};
	int val_t[Size*Size]={};
	laplacian(5, Size, Size, canny_1, lap, val_t);
	write_pgm_image("lap.pgm", Size, Size, lap);
	paint_histogram(Size, Size, lap, "histo_lap.pgm");

	unsigned char lap2[Size*Size]={};
	thresh_lap(Size, Size, lap, lap2, 207, 175);
	write_pgm_image("lap2.pgm", Size, Size, lap2);
	thresh_lap(Size, Size, lap, lap2, 206, 176);
	write_pgm_image("lap21.pgm", Size, Size, lap2);
	thresh_lap(Size, Size, lap, lap2, 205, 177);
	write_pgm_image("lap22.pgm", Size, Size, lap2);
	thresh_lap(Size, Size, lap, lap2, 204, 178);
	write_pgm_image("lap23.pgm", Size, Size, lap2);
	thresh_lap(Size, Size, lap, lap2, 200, 181);
	write_pgm_image("lap25.pgm", Size, Size, lap2);
	thresh_lap(Size, Size, lap, lap2, 198, 183);
	write_pgm_image("lap26.pgm", Size, Size, lap2);
	thresh_lap(Size, Size, lap, lap2, 202, 179);
	write_pgm_image("lap24.pgm", Size, Size, lap2);

	unsigned char lap3[Size*Size]={};
	zero_crossing(Size, Size, val_t, lap2, lap3);
	write_pgm_image("lap3.pgm", Size, Size, lap3);


	paint_histogram(Size, Size, lap2, "histo_lap2.pgm");

	exit(0);
	return 0;
}




	
