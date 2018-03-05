#ifndef __PWM_DRIVER_H__
#define __PWM_DRIVER_H__

#include "stm32f3xx_hal.h"
#include "stm32f3xx_hal_dma.h"
#include "stm32f3xx_hal_tim.h"
class pwm
{
	private:
		TIM_HandleTypeDef* m_timer;
	public:
		pwm(TIM_HandleTypeDef* timer);
		void set_timer_freq(uint32_t freq);
		void set_channel_duty_cycle(uint32_t Channel, 
								uint8_t duty_cycle);
		void start(uint32_t Channel);

};


#endif