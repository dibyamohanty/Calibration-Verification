//y rows x columns

#include <iostream>
#include <fstream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_cloud.h>
#include <cloud_viewer.h>
#include <pcl/impl/point_types.hpp>
#include <algorithm>
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/console/parse.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/passthrough.h>
#include <memory.h>
#include <string>

typedef pcl::PointXYZRGBA PointType;

using namespace pcl;
using namespace std;
using namespace Eigen;


pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGBA>);
pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_old (new pcl::PointCloud<pcl::PointXYZRGBA>);
//pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_filtered (new pcl::PointCloud<pcl::PointXYZRGBA>);
float x,y,z,v,*c_x,*c_y,*c_z,*n_x,*n_y,*n_z;

Matrix3Xf B_col; 

std::shared_ptr<pcl::visualization::PCLVisualizer> viewer2 (
    pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud);

int main (int argc, char** argv)
{
	/////////////////////////////////////Reading Data///////////////////////////////////

    pcl::PointCloud<pcl::PointXYZRGBA> cloud_in;

	ifstream in;
	ofstream out;
	in.open ("script.txt");

    // Fill in the cloud data
	float cloud_size;
    cloud_in.width    = 211;	//2464 211;
    cloud_in.height   = 90;	//2525 90;
	cloud_size = cloud_in.width * cloud_in.height;
    cloud_in.points.resize (cloud_in.width * cloud_in.height);
	uint8_t r,g,b;
    r = 0; 
	g = 0; 
	b = 0;
    unsigned int rgb = 0;
	int iPCDindex=0;
	
	while (in >> x)
	{
        in >> y >> z >> v;
		//in >> y >> v;
    	//cout << x << " , " << y << " , " << z << " , " << v << endl; 

		cloud_in.points[iPCDindex].x = x;
		cloud_in.points[iPCDindex].y = y;
		cloud_in.points[iPCDindex].z = z;
		r = v;                   
		g = v; 
		b = v;

		rgb = ((int)r) | ((int)g)  | ((int)b);
		//rgb = ((int)r) << 16 + ((int)g) << 8 + ((int)b);
		cloud_in.points[iPCDindex].rgba = rgb;

    	//cout << cloud_xyzrgb.points[iPCDindex].x << " , "<< cloud_xyzrgb.points[iPCDindex].y << " , "<<cloud_xyzrgb.points[iPCDindex].z << " , "<< endl;
		iPCDindex++;
	}

	cout << "saving PCD file .. " << endl;
	
	char file_name_old[100];
	sprintf(file_name_old, "b_old.pcd", 1);
	pcl::io::savePCDFileBinary (file_name_old, cloud_in );
	
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_old (new pcl::PointCloud<pcl::PointXYZRGBA>);
	
	if (pcl::io::loadPCDFile<pcl::PointXYZRGBA> ("b_old.pcd", *cloud_old) == -1) //* load the file
	{
		PCL_ERROR ("Couldn't read file b_old.pcd \n");
		system("PAUSE");

		return (-1);
	} 

	///////////////////////////Filtering///////////////////////////////////////////
	//Filtering
	pcl::PassThrough<pcl::PointXYZRGBA> pass;
	pass.setInputCloud (cloud);
	pass.setFilterFieldName ("x");
	pass.setFilterLimits (0.0, 0.2);
	//pass.setFilterLimitsNegative (true);
	pass.filter (*cloud);
	

	//////////////////////////////////////Processing//////////////////////////////////////

	//Regression
	float cloud_in_x=0,cloud_in_y=0,cloud_in_z=0,cloud_in_xx=0,cloud_in_xy=0,cloud_in_xz=0,cloud_in_yy=0,cloud_in_yz=0;
	for (int i=0; i<cloud_size; i++)
	{
		cloud_in_x=cloud_in_x+(cloud_in.points[i].x);
		cloud_in_y=cloud_in_y+(cloud_in.points[i].y);
		cloud_in_z=cloud_in_z+(cloud_in.points[i].z);
		cloud_in_xx=cloud_in_xx+((cloud_in.points[i].x)*(cloud_in.points[i].x));
		cloud_in_xy=cloud_in_xy+((cloud_in.points[i].x)*(cloud_in.points[i].y));
		cloud_in_xz=cloud_in_xz+((cloud_in.points[i].x)*(cloud_in.points[i].z));
		cloud_in_yy=cloud_in_yy+((cloud_in.points[i].y)*(cloud_in.points[i].y));
		cloud_in_yz=cloud_in_yz+((cloud_in.points[i].y)*(cloud_in.points[i].z));
	}

	Matrix3f A;
	Vector3f C;

	A << cloud_in_xx,cloud_in_xy,cloud_in_x, cloud_in_xy,cloud_in_yy,cloud_in_y, cloud_in_x,cloud_in_y,1;

	C << cloud_in_xz,  cloud_in_yz,  cloud_in_z;

	Vector3f B;
	B = A.inverse() * C;
	//cout<< B;
	c_x=B.row(0).data();
	c_y=B.row(1).data();
	c_z=B.row(2).data();
	
	//cout << *c_x << *c_y;
	Vector3f n_v(*c_x,*c_y,1), z_v(0,0,1);
	float dp, z_n, n_n, theta, rad;
	z_n = z_v.norm();
	n_n = n_v.norm();
	VectorXf cr, cr_sin2, cr_cos2, cr_sin, cr_cos;
	Matrix3f rot;
	cr = n_v.cross(z_v);
	dp = n_v.dot(z_v);
	theta = acos (dp/(z_n*n_n));
	rad=(3.14*theta)/180;
	//cout << theta << " " << rad << " ";

	rot << cos(rad),0,sin(rad), 0,1,0, -(sin(rad)),0,cos(rad);
	//cout << rot << " ";

	//3D to 2D rotated
	for (int i=0; i<cloud_size; i++)
	{
		Vector3f xyz, xyz_new;
		xyz << cloud_in.points[i].x,cloud_in.points[i].y,cloud_in.points[i].z;
		xyz_new = rot*xyz;
		cloud_in.points[i].x = xyz_new(0);
		cloud_in.points[i].y = xyz_new(1);
		cloud_in.points[i].z = xyz_new(2);
		//cout << cloud_in.points[i].x << " " << cloud_in.points[i].y << " " << cloud_in.points[i].z << " " ;
	}

	for (int i=0; i<cloud_size; i++)
	{
		cloud_in.points[i].z=0;
	}

	/*//Sorting
	for (int i=0; i<cloud_size; i++)
	{
		for (int j=i; j<cloud_size; j++)
		{
			float temp, tempx;
			temp = cloud_in.points[i].y;
			if (temp > cloud_in.points[j].y)
			{
				cloud_in.points[i].y = cloud_in.points[j].y;
				cloud_in.points[j].y = temp;
				tempx = cloud_in.points[i].x;
				cloud_in.points[i].x = cloud_in.points[j].x;
				cloud_in.points[j].x = temp;
			}
		}
	}
	*/

	float x_minn, y_minn, x_maxx, y_maxx;
	x_minn=cloud_in.points[0].x;
	y_minn=cloud_in.points[0].y;
	x_maxx=cloud_in.points[0].x;
	y_maxx=cloud_in.points[0].y;
	for (int i=0; i<cloud_size; i++)
	{
		if (cloud_in.points[i].x < (x_minn))
		{
			x_minn=cloud_in.points[i].x;
		}
			
		if (cloud_in.points[i].x > (x_maxx))
		{
			x_maxx=cloud_in.points[i].x;
		}

		if (cloud_in.points[i].y < (y_minn))
		{
			y_minn=cloud_in.points[i].y;
		}
		if (cloud_in.points[i].y > (y_maxx))
		{
			y_maxx=cloud_in.points[i].y;
		}
	}

	cout << x_minn << " " << y_minn << " " << x_maxx << " " << y_maxx << endl;

	for (int i=0; i<cloud_size; i++)
	{
		cloud_in.points[i].x=cloud_in.points[i].x + abs (x_minn);
		cloud_in.points[i].y=cloud_in.points[i].y + abs (y_minn);
	}

	//Search Window
	float x_max, x_min, y_max, y_min;
	x_min=cloud_in.points[0].x;
	y_min=cloud_in.points[0].y;
	x_max=cloud_in.points[0].x;
	y_max=cloud_in.points[0].y;
	for (int i=0; i<cloud_size; i++)
	{
		if (cloud_in.points[i].x >= (x_max))
		{
			x_max=cloud_in.points[i].x;
		}
		if (cloud_in.points[i].x <= (x_min))
		{
			x_min=cloud_in.points[i].x;
		}

		if (cloud_in.points[i].y >= (y_max))
		{
			y_max=cloud_in.points[i].y;
		}
		if (cloud_in.points[i].y <= (y_min))
		{
			y_min=cloud_in.points[i].y;
		}
	}

	cout << x_min << " " << y_min << " " << x_max << " " << y_max << endl;
	//cout << cloud_in.points[0].x << " " << cloud_in.points[0].y << " " << cloud_in.points[1].x << " " << cloud_in.points[1].y;

	ofstream datum("count.txt");
	float len, bre;
	int count_x=0, count_y=0, t_count=0;
	vector<float>windowcount_f;
	vector<vector<float>>windowcount;
	for (len=y_min-0.029; len <= y_max+0.029; len = len+0.002)
	{
		for (bre=x_min-0.029; bre <= x_max+0.029; bre = bre+0.002)
		{
			for (int i=0; i < cloud_size; i++)
			{
				if ((cloud_in.points[i].x > bre) && (cloud_in.points[i].x < bre+0.029) && (cloud_in.points[i].y > len) && (cloud_in.points[i].y < len+0.029))
				{
					t_count++;
				}
			}
			datum << t_count << "\t";
			//cout << t_count;
			windowcount_f.push_back(t_count);		
			//datum << windowcount_f.data();
			t_count=0;
			count_x++;
		}
		windowcount.push_back(windowcount_f);
		datum << endl;
		windowcount_f.erase(windowcount_f.begin(),windowcount_f.end());
		count_y++;
	}
	if ((count_x>0)||(count_y>0))
	{
		count_x = count_x/count_y;
	}
	count_y = count_y-1;
	cout << count_x << " " << count_y << " ";
	datum.close();

	//Derivatives and centre
	vector<float> minima_x;
	vector<float> minima_y;
	vector<float> minima_xc;
	vector<float> minima_yc;
	for (int i=1; i < count_y-1; i++) //rows
	{
		for (int j=1; j < count_x-1; j++) //columns
		{
			//cout << (windowcount[j+1][i] - windowcount[j][i]) << " " << (windowcount[j][i+1] - windowcount[j][i]) << " ";

			if (((windowcount[i+1][j] - windowcount[i][j]) > 0) && ((windowcount[i][j] - windowcount[i-1][j]) < 0) && ((windowcount[i][j+1] - windowcount[i][j]) > 0) && ((windowcount[i][j] - windowcount[i][j-1]) < 0))
			{
				//minima_x.push_back(y_min+(i*0.001));
				//minima_y.push_back(x_min+(j*0.001));
				minima_xc.push_back(((y_min - 0.029)+(i*0.002))+0.01499);
				minima_yc.push_back(((x_min - 0.029)+(j*0.002))+0.01499);
				//cout << (y_min+(i*0.02)) << " " << (x_min+(j*0.02)) << " " ;
				cout << (((x_min-0.029)+(j*0.002))+0.01499) << " " << (((y_min-0.029)+(i*0.002))+0.01499) << " " ;
			}
		}
	}

	



	/*//Sector search Window
	int p_count=0;
	float area = 0;
	vector<float>swindowcount_f;
	vector<vector<float>>swindowcount;
	for (int j=0; j < minima_xc.size(); j++)
	{
		for (theta=0; theta<360; theta=theta+5)
		{
			for (float rx = minima_x[j], float ry=minima_y[j]; rx < 0.029, ry < 0.029; rx=rx+0.001, ry=rx+0.001)
			{
				for (int i=0; i < cloud_size; i++)
				{
					if ((cloud_in.points[i].x > rx)&&(cloud_in.points[i].x > (rx*cos(theta)))&&(cloud_in.points[i].y > ry)&&(cloud_in.points[i].x > (ry*sin(theta))))
					{
						p_count++;
					}
				}
				area = 0.001
				swindowcount_f.push_back(p_count++);
				p_count=0;
			}
			swindowcount.push_back(swindowcount_f);
			swindowcount_f.erase(swindowcount_f.begin(),swindowcount_f.end());
		}
	}*/

	//Saving PCD file
	cout << "saving PCD file .. " << endl;
	
	char file_name[100];
	sprintf(file_name, "b.pcd", 1);
	pcl::io::savePCDFileBinary (file_name, cloud_in );
	
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGBA>);
	
	if (pcl::io::loadPCDFile<pcl::PointXYZRGBA> ("b.pcd", *cloud) == -1) //* load the file
	{
		PCL_ERROR ("Couldn't read file b.pcd \n");
		system("PAUSE");

		return (-1);
	} 
	

	if (pcl::io::loadPCDFile<pcl::PointXYZRGBA> ("b.pcd", *cloud) == -1) //* load the file
	{
		PCL_ERROR ("Couldn't read file b.pcd \n");
		system("PAUSE");

		return (-1);
	} 

	

	///////////////////////////////Visualize////////////////////////////////////
	//Calling visualizer
	std::shared_ptr<pcl::visualization::PCLVisualizer> vi;
	vi=viewer2(cloud);
	//vi=viewer2(cloud, n_v);
	while (!vi->wasStopped ())
	{
		vi->spinOnce (100);
		boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}
	system("PAUSE");

	
	return (0);
}


//Visualizer function
std::shared_ptr<pcl::visualization::PCLVisualizer> viewer2 (
   pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud)
{
  std::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
  viewer->setBackgroundColor (0, 0, 0);
  pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBA> rgb(cloud);
  viewer->addPointCloud<pcl::PointXYZRGBA> (cloud, rgb, "c.pcd");
  viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "c.pcd");
  return (viewer);
}



