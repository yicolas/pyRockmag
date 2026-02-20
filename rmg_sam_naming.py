"""
SAM file naming utilities with 8.3 format validation.
"""

def validate_sam_naming(site_id: str, sample_number: str, extension: str = "") -> tuple:
    """
    Validate SAM file naming according to 8.3 format.
    
    Parameters
    ----------
    site_id : str
        Site name (max 8 chars)
    sample_number : str
        Sample number/identifier
    extension : str
        File extension (max 3 chars after dot)
    
    Returns
    -------
    tuple : (is_valid, full_name, warning_message)
    """
    # Check site name length
    if len(site_id) > 8:
        return (False, "", f"Site name '{site_id}' exceeds 8 characters")
    
    # Build full filename
    if extension:
        # Remove leading dot if present
        ext = extension.lstrip('.')
        if len(ext) > 3:
            return (False, "", f"Extension '{extension}' exceeds 3 characters")
        full_name = f"{site_id}{sample_number}.{ext}"
    else:
        full_name = f"{site_id}{sample_number}"
    
    # Check total length (8.3 format)
    name_part, _, ext_part = full_name.partition('.')
    
    if len(name_part) > 8:
        return (False, full_name, 
                f"Filename '{name_part}' exceeds 8 characters (before extension)")
    
    if ext_part and len(ext_part) > 3:
        return (False, full_name,
                f"Extension '.{ext_part}' exceeds 3 characters")
    
    return (True, full_name, "")


def generate_sample_names(site_id: str, num_samples: int, 
                          use_integers: bool = True,
                          sequential: bool = True,
                          extension: str = "",
                          start_number: int = 1) -> list:
    """
    Generate sample names according to user preferences.
    
    Parameters
    ----------
    site_id : str
        Site identifier (max 8 chars)
    num_samples : int
        Number of samples to generate
    use_integers : bool
        If True, use integers (1, 2, 3). If False, use decimals (1.0, 2.0, 3.0)
    sequential : bool
        If True, sequential (1, 2, 3...). If False, sub-samples (1.1, 1.2, 1.3...)
    extension : str
        File extension (e.g., 'FRC', '0a')
    start_number : int
        Starting number for samples
    
    Returns
    -------
    list : List of (sample_name, is_valid, warning) tuples
    """
    samples = []
    
    for i in range(num_samples):
        if sequential:
            # Sequential: 1, 2, 3... or 1.0, 2.0, 3.0...
            number = start_number + i
            if use_integers:
                sample_num = str(number)
            else:
                sample_num = f"{number}.0"
        else:
            # Sub-samples: 1.1, 1.2, 1.3...
            sub_num = i + 1
            sample_num = f"{start_number}.{sub_num}"
        
        # Validate
        is_valid, full_name, warning = validate_sam_naming(site_id, sample_num, extension)
        samples.append((full_name, is_valid, warning))
    
    return samples


def get_sample_naming_interactive():
    """
    Interactive prompts for sample naming preferences.
    
    Returns
    -------
    dict : Dictionary with naming preferences
    """
    print("\n--- Sample Naming Configuration ---")
    
    # Number format
    print("\nNumber format:")
    print("  [1] Integers (1, 2, 3, ...)")
    print("  [2] Decimals (1.0, 2.0, 3.0, ...)")
    num_format = input("Choice [1]: ").strip() or "1"
    use_integers = (num_format == "1")
    
    # Sequential pattern
    print("\nNumbering pattern:")
    print("  [1] Sequential (1, 2, 3, ...)")
    print("  [2] Sub-samples (1.1, 1.2, 1.3, ...)")
    pattern = input("Choice [1]: ").strip() or "1"
    sequential = (pattern == "1")
    
    # Starting number
    try:
        start_num = int(input("\nStarting number [1]: ").strip() or "1")
    except ValueError:
        print("Invalid number, using 1")
        start_num = 1
    
    # Extension
    print("\nFile extension (max 3 characters):")
    print("  Examples: FRC (for FORC), 0a (standard), or blank for none")
    extension = input("Extension [blank]: ").strip()
    
    return {
        'use_integers': use_integers,
        'sequential': sequential,
        'start_number': start_num,
        'extension': extension
    }


def preview_sample_names(site_id: str, num_samples: int, naming_prefs: dict, 
                        preview_count: int = 5):
    """
    Show preview of sample names to user.
    
    Parameters
    ----------
    site_id : str
        Site identifier
    num_samples : int
        Total number of samples
    naming_prefs : dict
        Naming preferences from get_sample_naming_interactive()
    preview_count : int
        Number of examples to show
    """
    print(f"\n--- Sample Name Preview (first {min(preview_count, num_samples)}) ---")
    
    samples = generate_sample_names(
        site_id=site_id,
        num_samples=min(preview_count, num_samples),
        **naming_prefs
    )
    
    all_valid = True
    for full_name, is_valid, warning in samples:
        status = "✓" if is_valid else "✗"
        print(f"  {status} {full_name}")
        if not is_valid:
            print(f"      Warning: {warning}")
            all_valid = False
    
    if num_samples > preview_count:
        print(f"  ... ({num_samples - preview_count} more samples)")
    
    if not all_valid:
        print("\n⚠ WARNING: Some sample names exceed 8.3 format limits!")
        print("   Suggestion: Use shorter site name or different numbering")
        return False
    
    return True

