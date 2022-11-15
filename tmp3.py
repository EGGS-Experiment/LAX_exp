class MissingTrigger(Exception):
    pass

class ExternalTrigger(HasEnvironment):
    def build(self, trigger=None, t_timeout = 100*ms):
        self.setattr_device("core")
        self.trigger = trigger
        self.t_timeout = t_timeout

    def prepare(self):
        self.t_timeout_mu = self.core.seconds_to_mu(self.t_timeout)
        self.t_buffer_mu = self.core.seconds_to_mu(20*us)

    @kernel
    def wait_for_trigger(self):
        t_gate_open = now_mu()
        self.trigger._set_sensitivity(1)
        # Loop until all old (before current gate open) events are consumed, or
        # there is a timeout
        t_trig_mu = 0
        while True:
            # Wait for a trigger event for up to t_timeout_mu before returning
            t_trig_mu = rtio_input_timestamp(now_mu() + self.t_timeout_mu, self.trigger.channel)
            # If event if a timeout in the current gate period
            if t_trig_mu < 0 or t_trig_mu >= t_gate_open:
                break
        t_wall = self.core.get_rtio_counter_mu()
        at_mu(t_wall + self.t_buffer_mu)
        self.trigger._set_sensitivity(0)
        if t_trig_mu < 0:
            raise MissingTrigger()
        return t_trig_mu

class TTLin_block(EnvExperiment):
    def build(self):
        self.setattr_device("core")

        # Get all TTL out
        self.ttls = [ self.get_device("ttl"+str(i)) for i in range(4,64) ]
        self.start = self.get_device("ttl0")

    def prepare(self):
        self.lt = ExternalTrigger(self, self.start)
        self.lt.prepare()

    @kernel
    def run(self): # works, but 300 ns jitter
        self.core.reset()
        self.lt.wait_for_trigger()
        self.ttls[0].pulse(5*ms)