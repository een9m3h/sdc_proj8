/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <assert.h>     /* assert */

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	
	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

	// Set standard deviations for x, y, and theta
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];
	 

	// This line creates a normal (Gaussian) distribution for x, y and theta
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	
	for (int i = 0; i < num_particles; ++i) {
		//double sample_x, sample_y, sample_theta;
		
		// TODO: Sample  and from these normal distrubtions like this: 
		//	 sample_x = dist_x(gen);
		//	 where "gen" is the random engine initialized earlier.
		
		Particle p = Particle();
		p.id = i;
		p.x	= dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		
		particles.push_back(p);
		
		/* sample_x = dist_x(gen);
		 sample_y = dist_y(gen);
		 sample_theta = dist_theta(gen);	 
		 
		 // Print your samples to the terminal.
		 cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " << sample_theta << endl;
		 */
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//assert(velocity > 0.0);
	
	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

	// Set standard deviations for x, y, and theta
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];
	
	
	for (int i = 0; i < num_particles; ++i) {
		Particle p = particles[i];
		
		// This line creates a normal (Gaussian) distribution for x, y and theta
		normal_distribution<double> dist_x(p.x, std_x);
		normal_distribution<double> dist_y(p.y, std_y);
		normal_distribution<double> dist_theta(p.theta, std_theta);
		
		p.theta = dist_theta(gen);
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		
		//avoid division by zero
		if (fabs(yaw_rate) > 0.001) {
			p.x = p.x + velocity/yaw_rate * ( sin (p.theta + yaw_rate*delta_t) - sin(p.theta));
			p.y = p.y + velocity/yaw_rate * ( cos(p.theta) - cos(p.theta + yaw_rate*delta_t) );
		}
		else {
			p.x = p.x + velocity*delta_t*cos(p.theta);
			p.y = p.y + velocity*delta_t*sin(p.theta);
		}
		
		particles[i] = p;
	}

}

std::vector<LandmarkObs> ParticleFilter::dataAssociation(Particle &p, double sensor_range, std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	std::vector<LandmarkObs> selPredicted;
	std::vector<double>  sense_x, sense_y;
	std::vector<int> associations;
	
	for(int i = 0; i < observations.size(); i++){
		double min_dist = sensor_range;
		LandmarkObs currPredicted;
		for(int j = 0; j < predicted.size(); j++){
			double curr_dist = calcEuclidDist(observations[i], predicted[j]);
			if( curr_dist < min_dist){
				min_dist = curr_dist;
				observations[i].id = predicted[j].id;
				currPredicted = predicted[j];
			}			
		}
		selPredicted.push_back(currPredicted);
		sense_x.push_back(observations[i].x);
		sense_y.push_back(observations[i].y);
		associations.push_back(observations[i].id);
	}
	
	SetAssociations(p, associations, sense_x, sense_y);
	
	return selPredicted;	
}

double ParticleFilter::calcParticleWeight(double std_landmark[], std::vector<LandmarkObs> a, std::vector<LandmarkObs> b){
	
	assert (a.size() == b.size());
	double totalWeight = 0.0;
	
	for(int i = 0; i < a.size(); i++){
		totalWeight *= calcWeight(std_landmark, a[i], b[i]);
	}
	
	return totalWeight;
}

double ParticleFilter::calcWeight(double std_landmark[], LandmarkObs a, LandmarkObs b){
	
	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];
	
	
	// calculate normalization term
	double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));

	// calculate exponent
	double exponent= (pow((a.x - b.x),2))/(2 * pow(sig_x,2)) + (pow((a.y - b.y),2))/(2 * pow(sig_y,2));

	// calculate weight using normalization terms and exponent
	double weight= gauss_norm * exp(-exponent);
	
	return weight;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for(int i = 0; i < num_particles; i++){
		Particle p = particles[i];
		
		//make car observation copy to be translated into map co-ordinates by particle position
		std::vector<LandmarkObs> localObs(observations.begin(), observations.end());
		
		// Translate observation co-ordinates from car to map based on particle orientation and position
		mapTransate(p, localObs);
		
		//filter map features based on distance to particle and sensor range
		std::vector<LandmarkObs> predicted = rangeFilter(sensor_range, p, map_landmarks);
		
		//make associations
		std::vector<LandmarkObs> assocLandmarks = dataAssociation(p, sensor_range, predicted, localObs);
		
		p.weight = calcParticleWeight(std_landmark, assocLandmarks, localObs);
		
	}
	
	//normalize weights to probabilities
	double totalWeights = 0;
	for(int i = 0; i < num_particles; i++){
		totalWeights += particles[i].weight;
	}
	for(int i = 0; i < num_particles; i++){
		particles[i].weight /= totalWeights;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,num_particles-1);
	std::vector<Particle> newParticles;
	
	double max_weight = 0;
	for(int i = 0; i < num_particles; i++){
		if(particles[i].weight > max_weight)
			max_weight = particles[i].weight;
	}
	
	int index = distribution(generator);
	
	double beta = 0.0;
	for(int i = 0; i < num_particles; i++){
		beta += (distribution(generator)/num_particles)*2.0*max_weight;
		while(beta > particles[index].weight){
			beta -= particles[index].weight;
			index = (index + 1 ) % num_particles;
		}
		newParticles.push_back(particles[index]);
	}
	particles = newParticles;
	

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

double ParticleFilter::calcEuclidDist(LandmarkObs a, LandmarkObs b)
{
	double x = a.x - b.x; //calculating number to square in next step
	double y = a.y - b.y;
	double dist;

	dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
	dist = sqrt(dist);                  

	return dist;
}

double ParticleFilter::calcEuclidDist(LandmarkObs a, Map::single_landmark_s b)
{
	double x = a.x - b.x_f; //calculating number to square in next step
	double y = a.y - b.y_f;
	double dist;

	dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
	dist = sqrt(dist);                  

	return dist;
}

void ParticleFilter::mapTransate(const Particle& p, std::vector<LandmarkObs> &observations){
	
	
	for(int i = 0; i < observations.size(); i++){
		
		double x_map, y_map;
		
		LandmarkObs obs = observations[i];
	
		///* transform to map x coordinate
		x_map= p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);

		///* transform to map y coordinate
		y_map= p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);
		
		observations[i].x = x_map;
		observations[i].y = y_map;
	}
}

std::vector<LandmarkObs> ParticleFilter::rangeFilter(double sensor_range, Particle &p, const Map &map_landmarks){
	
	std::vector<LandmarkObs> retainedLandmarks;
	
	LandmarkObs particlePosition = {0, p.x, p.y};
	
	for (int i = 0 ; i < map_landmarks.landmark_list.size(); i++){
		LandmarkObs newLandmark;
		if(calcEuclidDist(particlePosition, map_landmarks.landmark_list[i]) <= sensor_range){
			newLandmark.id = map_landmarks.landmark_list[i].id_i;
			newLandmark.x = map_landmarks.landmark_list[i].x_f;
			newLandmark.y = map_landmarks.landmark_list[i].y_f;
			retainedLandmarks.push_back(newLandmark);
		}
			 
	}
	return retainedLandmarks;
	
}