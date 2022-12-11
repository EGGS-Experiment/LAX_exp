"""
DMA Decorator

Combines multiple kernel sequences into a single sequence and records it onto DMA.
# todo: finish
"""
def autoreload(func):
    """
    A decorator for non-kernel functions that attempts to reset
    the connection to artiq_master if we lost it.
    """

    def inner(self, *args, **kwargs):

        # run decorated function
        try:
            res = func(self, *args, **kwargs)
            return res

        # handle connection errors
        except (ConnectionAbortedError, ConnectionResetError) as e:

            # attempt to reconnect to master and reset client objects
            try:
                # reset connection
                print('Connection aborted, resetting connection to artiq_master...')
                self.reset_connection()

                # retry decorated function
                res = func(self, *args, **kwargs)
                return res
            except Exception as e:
                raise e

        # handle RTIOUnderflows
        except RTIOUnderflow as e:
            # add slack
            self.core.break_realtime()

            # retry decorated function
            res = func(self, *args, **kwargs)
            return res

        # handle RTIOOverflows
        except RTIOOverflow as e:
            # reset core & add more slack
            self.core.reset()
            self.core.break_realtime()

            # retry decorated function
            res = func(self, *args, **kwargs)
            return res

    return inner
