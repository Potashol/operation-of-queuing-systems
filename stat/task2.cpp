#include<iostream>
#include <cmath>
#include <vector>
#include <random>
#include <functional>
namespace stat_mod_3 {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution <> urd(0, 1);
	double uniform()
	{
		return urd(gen);
	}
	double exp_random_value(double lambda)
	{
		return(-1.0 / lambda) * log(uniform());
	}
	
	// Построение потока входных заявок
	std::vector<double>get_request_times(double modeling_time, double lambda)
	{
		std::vector<double>request_times;
		double time  = exp_random_value(lambda);
		do
		{
			request_times.push_back(time);
			time+= exp_random_value(lambda);
		}while(time < modeling_time);
		return request_times;
	}

	//Имитирование работы СМО:
 
		//Модель каналов обслуживания

	class service_channels
	{
	private:
		std::vector <double> channels_service_end_times;
		int channels_count;
	public:
		service_channels(int count)
		{
			for(int channel = 0; channel < count; channel++)
			{
				channels_service_end_times.push_back(0);
			}
			channels_count=count;
		}
		double get_service_end_time(int channel_index)
		{
			return channels_service_end_times[channel_index];
		}
		int get_closest_service_end_time_channel()
		{
			double closest_time = INFINITY;
			int channel_index = 0;
			for(int i = 0; i <channels_count; i++)
			{
				if (channels_service_end_times [i] < closest_time)
				{
					closest_time = channels_service_end_times[i];
					channel_index = i;
				}
			}
		return channel_index;
		}
		void update_channel_service_end_time (int channel_index, double new_time)
		{
			channels_service_end_times[channel_index] = new_time;
		}
	};

	struct smo_event
	{
		double time;
		int state_change_factor;//+1 or −1
		smo_event(double time, int state_change_factor)
		:time(time), state_change_factor(state_change_factor){}
		};
		class smo_events_aggregator
		{
		private:
			std::vector<smo_event*>events;
			void add_event(double time, int state_change_factor)
			{
				smo_event*e=new smo_event(time,state_change_factor);
				int size=events.size();
				int insert_index=size;
				for(int i = 0; i<size; i++)
				{
					if(events[i] -> time > time)
					{
						insert_index=i;
						break;
					}
				}
				events.insert(events.begin() + insert_index,e);
			}
		public:
		void add_request (double begin_time, double end_time)
		{
			add_event(begin_time,1);
			add_event(end_time,-1);
		}
		int get_state_at_time(double time)
		{
			int state=0;
			for(smo_event*e:events)
			{
				if(e -> time > time)
				{
					break;
				}
				state+= e-> state_change_factor;
			}
			return state;
		}
		double get_state_total_time(int state)
		{
			int count=0;
			double prev_event_time=0;
			double state_total_time=0;
			for(smo_event*e:events)
			{
				if(count==state)
				{
					state_total_time+=(e -> time-prev_event_time);
				}
				count+= e-> state_change_factor;
				prev_event_time= e-> time;
			}
			return state_total_time;
		}
		~smo_events_aggregator()
		{
			for(int i=0; i < events.size(); i++)
			{
				delete events[i];
			}
			events.clear();
		}
	};

	//Моделирование процесса обслуживания

	double smo_modeling_process(int channels_count, std::vector <double>request_times, std::function <double()> get_service_time, int max_queue_length, double modeling_time, int smo_state)
	{
		service_channels channels (channels_count);
		smo_events_aggregator queue;
		smo_events_aggregator smo;
		int count = 0; 
		for(double request_time :request_times)
		{
			double request_exit_time = request_time;
			int channel = channels.get_closest_service_end_time_channel();
			double closest_service_end_time = channels.get_service_end_time(channel);
			if(closest_service_end_time < request_time)
			{
				//go to service immediately
				double service_time = get_service_time();
				request_exit_time = closest_service_end_time + get_service_time();
				channels.update_channel_service_end_time (channel,request_time + get_service_time());
			}
			else
			{
				//go to queue then to service
				++count; 
				queue.add_request(request_time, closest_service_end_time);
				double service_time = get_service_time();
				channels.update_channel_service_end_time(channel, closest_service_end_time + service_time);
				request_exit_time = closest_service_end_time + service_time;
			}
			//skip rejected requests
			if(request_time < request_exit_time)
			{
				smo.add_request(request_time, request_exit_time);
			}
		}
		//std::cout << max_chanels_count<< std::endl;
		return double(count) / request_times.size();
	}

	const int MIN_ITER_COUNT = 50;
	const int n = 3;
	const double N = 7;
	double expected_value(std::function<double()>rand_distr,double precision, int min_iter_count = MIN_ITER_COUNT)
	{
		double sum = 0;
		double square_sum = 0;
		double n = 0;
		bool is_precision_obtained = false;
		while(!is_precision_obtained)
		{
			double r = rand_distr();
			sum+= r;
			square_sum+=(r*r);
			n++;
			if(n < MIN_ITER_COUNT)
				continue;
			double dispersion = (1.0/(n-1))*square_sum-(1.0/(n * (n-1))) * sum * sum;
			double estimated_iter_count = 9 * dispersion / (precision * precision);
			is_precision_obtained = (n > estimated_iter_count);
		}
		return sum/n;
	}


	double operating_characteristic_numerical(int lambda, int m)
	{
		double mu = 60.0 / m;
		auto get_service_time = [&]() -> double
		{
		return exp_random_value(mu);
		};
		auto get_smo_state_prob = [&]() -> double
		{
		double modeling_time = 1000;
		auto request_times=get_request_times(modeling_time, lambda);
		return smo_modeling_process(n, request_times , get_service_time, N-n, modeling_time , 1);
		};
		return expected_value(get_smo_state_prob,0.01);
	}


	int fact(int n)
	{
		int sum = 1; 
		for (int i = 1; i <= n; ++i)
			sum *= i; 
		return sum; 
	}
	double operating_characteristic_analytical(int lambda, int m)
	{
		double mu = 60.0 / m;
		double ro = lambda / mu;
		double ro_n = pow(ro, n); 
		double p0 = 0.0;
		double sum = 0; 
		for (int k = 0; k < n; ++k) {
			p0 += pow(ro, k) /fact(k) + (ro_n/fact(n))*(1/(1 - ro/n));
		}
		return ((p0 * pow(ro, n + 1)) / (fact(n) * n)) * (1 / (pow((1 - ro / n), 2)));
	}
}

int main()
{
	std::cout << "Numerical:" << stat_mod_3::operating_characteristic_numerical(4, 10) << std::endl;
	std::cout << "Analytical:" << stat_mod_3::operating_characteristic_analytical(4, 10) << std::endl << std::endl;
}