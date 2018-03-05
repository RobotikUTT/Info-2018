#include "pwm.h"

pwm::pwm(TIM_HandleTypeDef* timer)
{
  m_timer = timer;
}

void pwm::set_timer_freq(uint32_t freq)
{
  uint32_t sys_freq = HAL_RCC_GetPCLK1Freq();

  m_timer->Init.Prescaler = 0;
  m_timer->Init.CounterMode = TIM_COUNTERMODE_UP;
  m_timer->Init.Period = (uint16_t)(2*sys_freq/(freq));
  m_timer->Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  m_timer->Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_PWM_Init(m_timer) != HAL_OK)
  {
    _Error_Handler(__FILE__, __LINE__);
  }
}

void pwm::set_channel_duty_cycle(uint32_t Channel, 
                                 uint8_t duty_cycle)
{
	TIM_OC_InitTypeDef sConfigOC;
	sConfigOC.OCMode = TIM_OCMODE_PWM1;
	sConfigOC.Pulse = m_timer->Init.Period*duty_cycle/0xFF;
	sConfigOC.OCPolarity = TIM_OCPOLARITY_HIGH;
	sConfigOC.OCFastMode = TIM_OCFAST_DISABLE;
	if (HAL_TIM_PWM_ConfigChannel(m_timer, &sConfigOC, Channel) != HAL_OK)
	{
	  _Error_Handler(__FILE__, __LINE__);
	}
  start(Channel);
}

void pwm::start(uint32_t Channel)
{
  HAL_TIM_PWM_Start(m_timer,Channel);
}