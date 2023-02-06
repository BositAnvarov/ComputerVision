// Authors: Bosit Anvarov, Hanat Samatar, and Kush Chopra
// Topic: Road Sign Detection (further details in the writeup)

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>
#include "opencv2/imgcodecs.hpp"
#include <iostream>
using namespace cv;
using namespace std;


// enum to label the most important colors we are trying to find in the image
// Traffic Sign Color Correspondance:
//			RED - any sign that must be followed. For example: Stop Sign, Yield Sign, Do Not Enter, etc...
//			YELLOW - any sign that is a warning. For example: Construction Sign, Merge Sign, Traffic Light Ahead, etc..
//			FLUORSCENT - similar to warning, but with high emphasis. For example: Road signs in school driveway, bicyclist, etc..
//			BLUE - this sign represents there is a Hospital near
enum Color { RED, YELLOW, FLUORESCENT, BLUE };

// These are global variables to store the color within a given range.
// Since the most common color is not always specifc, giving them a upper and lower 
// limits will give the approximate color.
// upper and lower bounds for red signs
const Scalar red_low = Scalar(0, 100, 100);
const Scalar red_high = Scalar(10, 255, 255);

// upper and lower bounds for yellow signs
const Scalar yellow_low = Scalar(15, 100, 100);
const Scalar yellow_high = Scalar(34, 255, 255);

// upper and lower bounds for fluorscent signs
const Scalar fluorscent_low = Scalar(35, 100, 100);
const Scalar fluorscent_high = Scalar(50, 255, 255);

// upper and lower bounds for blue signs
const Scalar blue_low = Scalar(95, 100, 100);
const Scalar blue_high = Scalar(120, 255, 255);


// basic Sign struct to initialize the color, 
// rectangle, and the points in the rectangle
struct Sign
{
	Color color;
	Point position;
	Rect rect;

	Sign(Color color, Rect rect, int x, int y)
	{
		this->color = color;
		this->rect = rect;
		this->position = Point(x, y);
	}

};

// signs list to keep track of the number of different colors found in the image
vector<Sign> signs;


// largestContours() - finds and returns the biggest contours within the image
// precondition: the contours are fixed and given in the parameter
// postconditions: according the the maximum area generated by each contours
//				   we find the largest one and store it inside our result
Rect largestContours(const double& area, const vector<vector<Point>>& conts)
{
	// keep track of the maximums 
	Rect res;
	double maxArea = area;

	// iterate over each contour
	for (int i = 0; i < conts.size(); i++)
	{
		// store the area of the current contour
		double area = contourArea(conts[i]);

		// compare current contourArea with old maxContour area
		if (area > maxArea)
		{
			// assign their new values
			maxArea = area;
			res = boundingRect(conts[i]);
		}
	}

	return res;
}


// getSigns() - gets the largest contour for specific given color and places it inside Signs vector
// precondition: the input image exists in the directory
// postconditions: the input image gets converted to hsv colospace,
//				   using the inRange() function, it turns the background of 
//				   the image into black while turning the color we are looking 
//				   for into white, after that the contour can be detected
void getSigns(const Mat& input, Scalar low, Scalar high, Color color)
{
	// create matricies to store new data
	Mat mask, img;

	// convert rgb image to hsv
	// mask the hsv image using the given upper and lower bounds
	cvtColor(input, img, COLOR_BGR2HSV);
	inRange(img, low, high, mask);

	// 2d vector to keep track of contour points 
	vector<vector<Point>> contours;

	// find the contour points from the mask and record it
	findContours(mask, contours, RETR_EXTERNAL, CHAIN_APPROX_NONE);

	// get the biggest contours 
	Rect rect = largestContours(0, contours);

	// specify the size of the rectangle so it 
	// only detects the ones that are large enough
	if (rect.width > 200 && rect.height > 200)
	{
		// store all the biggest contours in the vector
		signs.emplace_back(color, rect, rect.x + rect.width / 2, rect.y + rect.height / 2);
	}

}

// showSigns() - this checks which of the signs are residing inside the image and returns the result
// precondition: input image is the original image
// postconditions: using the given vector length, we iterate through each Sign in the vector
//				   and draw the appropriate rectangle along with its name for each case that occurs
Mat showSigns(const Mat& input, const int& length)
{
	Mat res = input.clone();
	for (int i = 0; i < length; i++)
	{
		// RED, YELLOW, GREEN, FLUORESCENT, BLUE
		// draw the appropriate rectangle for each color
		// write the appropriate name for each color
		switch (signs[i].color) {

		case RED:
			rectangle(res, signs[i].rect.tl(), signs[i].rect.br(), Scalar(0, 255, 0), 5);
			putText(res, "REGULATORY", Point(signs[i].rect.x, signs[i].rect.y), FONT_HERSHEY_TRIPLEX, 1, Scalar(0, 0, 255), 2);
			break;

		case YELLOW:
			rectangle(res, signs[i].rect.tl(), signs[i].rect.br(), Scalar(255, 255, 0), 5);
			putText(res, "WARNING", Point(signs[i].rect.x, signs[i].rect.y), FONT_HERSHEY_TRIPLEX, 1, Scalar(0, 255, 255), 2);
			break;

		case FLUORESCENT:
			rectangle(res, signs[i].rect.tl(), signs[i].rect.br(), Scalar(255, 255, 0), 5);
			putText(res, "HIGH EMPHASIS WARNING", Point(signs[i].rect.x, signs[i].rect.y), FONT_HERSHEY_TRIPLEX, 1, Scalar(100, 100, 0), 2);
			break;

		case BLUE:
			rectangle(res, signs[i].rect.tl(), signs[i].rect.br(), Scalar(0, 0, 255), 5);
			putText(res, "GUIDANCE", Point(signs[i].rect.x, signs[i].rect.y), FONT_HERSHEY_TRIPLEX, 1, Scalar(255, 0, 0), 2);
			break;

		default:
			break;
		}
	}

	return res;
}


// Main 
// replace the imread() content with a different image in order to run a differnt image
// make sure the image is inside the directory of this program
int main(int argc, char* argv[])
{
	// read the given image
	Mat res = imread("test5.jpg");

	// check if the colors in the enum are present in the given input image
	// if so, put them in the signs vector
	getSigns(res.clone(), red_low, red_high, RED);
	getSigns(res.clone(), yellow_low, yellow_high, YELLOW);
	getSigns(res.clone(), blue_low, blue_high, BLUE);
	getSigns(res.clone(), fluorscent_low, fluorscent_high, FLUORESCENT);

	// get the resulting image which has the road signs detected if there is a match
	Mat output = showSigns(res, signs.size());

	// show the output on the screen
	// write the output in the directory
	imshow("output", output);
	waitKey(0);
	imwrite("output.jpg", output);

	return 0;
}