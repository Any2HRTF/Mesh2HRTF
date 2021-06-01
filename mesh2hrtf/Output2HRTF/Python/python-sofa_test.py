import sofa
print (sofa)

print(sofa.conventions.implemented())

HRIR_path = "free_field_HRIR.sofa"
measurements = 5
data_length = 1000
max_string_length = 128

# we will add this one later
receivers = 2


HRIR = sofa.Database.create(HRIR_path, "SimpleFreeFieldHRIR",
                            dimensions={"M": measurements, "N": data_length, "S": max_string_length})


try: HRIR.Listener.initialize()
except Exception as e: print("Initializing without Position error: {0}".format(e))

try: HRIR.Listener.initialize(fixed=["Position"])
except Exception as e: print("Initializing with incomplete convention requirements: {0}".format(e))

HRIR.Listener.initialize(fixed=["Position", "View", "Up"])
print("Successfully initialized Listener with fixed Position, View and Up.")


HRIR.Source.initialize(variances=["Position"])
HRIR.Source.initialize_coordinates(variances=["View"])

HRIR.Source.initialize_coordinates(fixed=["Position"]) # no change, message or output.


HRIR.Receiver.initialize(fixed=["Position"], count=receivers)

try: HRIR.Emitter.initialize(fixed=["Position"], count=2)
except Exception as e: print("Initializing with Emitter count not allowed in convention error: {0}".format(e))
HRIR.Emitter.initialize(fixed=["Position"])


print(HRIR.Data.Type)
HRIR.Data.initialize(variances=["Delay"])


HRIR.Room.Type = "shoebox"
HRIR.Room.initialize(variances=["CornerA", "CornerB"])

HRIR.Room.create_attribute("Location", "various recording locations")
HRIR.Room.create_variable("Temperature", ("M",))
HRIR.Room.Temperature.Units = "kelvin"
HRIR.Room.Temperature = 150
HRIR.Room.create_string_array("Description", ("M", "S"))

print(HRIR.Room.Location)
print(HRIR.Room.Temperature.get_values(), HRIR.Room.Temperature.Units)


HRIR.Metadata.set_attribute('GLOBAL_ApplicationName', 'Mesh2HRTF')

print("Attributes and metadata")
HRIR.Metadata.dump()


print("Dimensions")
HRIR.Dimensions.dump()


print("Variables")
HRIR.Variables.dump()


HRIR.close()